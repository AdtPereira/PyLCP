# Plano — Pipeline híbrido: parâmetros internos (FEM/COMSOL) + retorno pelo solo (analítico)

---

## STATUS — implementado para a Configuração 2

| Entregue | Arquivo |
|---|---|
| **Parte A** — retorno pela terra vê o raio externo do **tubo** (`_cable_external_geometry`): quando há `enclosure`, termo próprio de `Zg`/`Pg` usa `D2/2` e a profundidade do eixo do tubo | `mtl_main/strategy.py` (`SingleCoreCableInHDPEStrategy` + `...WithECCInHDPEStrategy`) |
| **Parte B** — `InternalParametersFromFEM`: `Zi` (interp. log-log) + `C` nodal → dict de 6 chaves, tilado `kron(I₃,·)` | `analytical_forms/single_core_cable.py` |
| Leitor `get_scc_internal_capacitance_matrix_combined()` (matriz nodal `C` crua) | `utils/comsol_data.py` |
| Orquestrador `build_pul_matrices(internal_source='analytical'\|'fem', ...)` | `analytical_forms/single_core_cable.py` |
| `andreata_case2`: cenário `'fem'` (Zi/Yi do COMSOL + terra analítica c/ raio do tubo), decomposição modal sobre ele, validação vs. MATLAB, 4 gráficos de comparação + Figs 5.8–5.10 | `testData/andreata_case2/` |

**Validação `andreata_case2` — cenário `'fem'` vs. MATLAB (FEM de Andreata):**

| grandeza | erro rel. máx. (norma ∞) |
|---|---|
| `Z'` (série completa) | **0,40 %** |
| `Y'` (shunt completa) | **1,85 %** |
| `Zg` (retorno pela terra) | **0,07 %** |
| `Pg` (coef. de potencial) | **0,81 %** |

O caminho FEM-híbrido reproduz a referência FEM de Andreata a **< 2 %** — muito
acima do GMD de Lafaia (mantido nos cenários `1/2/3` só para comparação).
Modal: rótulos limpos, `classif = 1,000`, `PASSIVE`; a oscilação de AF dos modos
coaxiais (artefato de degenerescência do GMD) **desaparece** com o `Zi` FEM.

Regressão: `andreata_case1/3/4` inalterados (case4 usa `scc-flat-ecc`, não HDPE).

**Pendências:** Configs 4–5 (bloco heterogêneo 3×3 do tubo compartilhado — a
Parte A já cobre a geometria; falta o `Zi` FEM 3×3); solo dependente da
frequência; pequeno offset de ~2 % em `C₂₂` (diferença entre o FEM do COMSOL
e o FEM próprio de Andreata, não do pipeline).

---

## 1. Motivação e diagnóstico

### Formulação de 2 termos (Ametani 2015) — a do pyLCP

O pyLCP decompõe os parâmetros PUL de um cabo enterrado em **dois** termos:

```
Z'(f) = Zi(f) + Zg(f)          Y'(f) = ( Yi(f)^-1 + Yg(f)^-1 )^-1
        \_ interno _/  \_ solo _/
```

- **`Zi` / `Yi`** — tudo que **não** é retorno pela terra: impedância dos
  condutores (pelicular/proximidade) **e** a indutância/capacitância externa da
  isolação da configuração até a sua superfície externa. Para Config. 1 essa
  superfície é a cobertura do SCC (`r4`); para Config. 2 é a **superfície
  externa do tubo HDPE** (`D2/2`), com ar + HDPE dentro do termo interno.
- **`Zg` / `Yg`** — retorno pela terra a partir dessa superfície, por expressões
  integrais (Magalhães/Xue) ou fechadas (De Conti/Duarte/Alipio 2023), que
  **só enxergam a posição e a maior dimensão externa** de cada objeto
  (Andreata eq. 4.27: *"x_jk igual ao raio externo dos cabos ou tubos"*).

> Andreata (§4.1, §4.4) escreve `Z' = Zi + Ze + Zg`, quebrando o termo interno
> em condutores (`Zi`, Bessel) + isolação externa (`Ze`, FEM). É a **mesma
> decomposição**: o `Zi` de 2 termos do pyLCP = `Zi + Ze` de Andreata. O FEM do
> COMSOL entrega justamente esse `Zi + Ze` já somado (ver §5).

**Premissa confirmada:** a fronteira externa do modelo COMSOL da Config. 2 é a
superfície externa do tubo HDPE (não há caixa de ar). Logo `Zi_fem` e `Zg`
casam exatamente no círculo `D2/2`, sem sobreposição nem lacuna.

### O que o pyLCP faz hoje

| | Config. 1 (`andreata_case1`) | Config. 2 (`andreata_case2`) |
|---|---|---|
| **interno** | analítico (Bessel núcleo/blindagem + forma fechada) — correto | analítico **+ truque GMD Case 3.1 de Lafaia**: o duto é dobrado numa isolação equivalente de permissividade `ε_a` na blindagem (`_override_sheath_insulation`) |
| **solo** — raio externo do termo próprio | `0.024 m` = raio externo do SCC (`r4`) — **correto** | `0.024 m` (cenário `3`/GMD) ou `r7` (cenário `2`/ERS) — **o tubo HDPE (`D2/2 = 0.0578 m`) é ignorado no cenário GMD** |

### Problemas

1. **Retorno pelo solo com raio errado (Parte A).**
   `SingleCoreCableInHDPEStrategy._cable_distance_matrices` seleciona o
   componente de maior raio externo por `data['radius'][1] + insulation.thickness`
   — **nunca olha `data['enclosure']`**. O tubo HDPE some do termo próprio de
   `d_matrix`/`D_matrix`. Verificado:
   `flat_hdpe_enclosed_model` → `diag(d_matrix_ground_return) = [0.024, 0.024, 0.024]`
   enquanto `enclosure['outer_radius'] = 0.0578`.

2. **Interno via GMD Case 3.1 (Parte B).**
   O truque de permissividade equivalente é uma aproximação (Lafaia 2015). O
   pyLCP **já lê** o `Zi`/`Yi` FEM do COMSOL para a geometria excêntrica exata
   da Config. 2 — hoje só como *overlay* de comparação, não como fonte de `Z'/Y'`.

### Objetivo

Reestruturar o pipeline para **separar as duas fontes**:

- **interno**: *pluggável* — `'analytical'` (geometrias canônicas: SCC nu,
  SCC em tubo metálico) **ou** `'fem'` (COMSOL, para geometrias não canônicas:
  SCC em tubo HDPE — Config. 2 — e adiante Configs 3–5);
- **solo**: **sempre analítico** (`magalhaes_xue` / `deconti`), usando a
  **maior dimensão externa real** da configuração (SCC → `r4`; tubo → `D2/2`).

Resultado: Config. 2 sem o GMD de Lafaia.

---

## 2. Princípio de projeto

O ponto de composição **já existe e não muda**:
`PerUnitParameters.quasi_tem_approx_matrices(internal_matrices, earth_return)`
faz `Zs = Zi + Zg` e `Ysh = jω (Pi + Pg)^-1`, lendo `Zi`, `Pi` e `block_sizes`
de um **dict `internal_matrices`**. Basta:

- **(A)** consertar o raio externo que `earth_return_parameters` usa (via
  `d_matrix_ground_return` etc.);
- **(B)** produzir o dict `internal_matrices` a partir do FEM, com a **mesma
  forma** que `InternalPerUnitParameters.matrices()` devolve hoje;
- **(C)** um orquestrador fino que escolhe a fonte interna e chama a composição.

Nenhuma mudança em `quasi_tem_approx_matrices`, `ModalDecomposition`,
`ModalPropagationPlotter`, `check_pul_passivity`.

---

## 3. Parte A — raio externo correto no retorno pelo solo

### 3.1 Onde

`mtl_main/strategy.py`:
- `SingleCoreCableInHDPEStrategy._cable_distance_matrices` (~linha 586)
- `SingleCoreCableWithECCInHDPEStrategy._cable_distance_matrices` (mesma cópia)

`SingleCoreCableStrategy` (SCC nu, ~linha 360) **não muda** — lá não há `enclosure`.

### 3.2 Mudança

Ao agrupar condutores por `center_point` e escolher o representante do cabo,
considerar o `enclosure` (tubo):

```python
def get_outer_radius(conductor_tuple):
    data = conductor_tuple[1]
    encl = data.get('enclosure')
    if encl and encl.get('outer_radius'):
        return encl['outer_radius']                    # tubo domina
    insul = (data.get('insulation') or {}).get('thickness', 0)
    return data['radius'][1] + insul
```

e no termo próprio (`n_tag == m_tag`):

```python
encl = n_cable.get('enclosure')
s = (encl['outer_radius'] if encl and encl.get('outer_radius')
     else n_cable['radius'][1] + (n_cable.get('insulation') or {}).get('thickness', 0))
```

**Profundidade do termo próprio**: quando o tubo é o objeto representativo, usar
o **centro do tubo** (`enclosure['center_point'][1]`, ex. `-1.1732`) em vez do
centro do cabo (`-1.2`) nos termos `images_vertical_distance_matrix` e `D_matrix`
próprios. Os termos **mútuos** (centro-a-centro, 0,2 m) **não mudam** — os tubos
não se tocam (`0.2 > 2 · 0.0578`).

### 3.3 Alternativa mais limpa (recomendada)

Em vez de espalhar a lógica pelo laço, o strategy calcula **um campo por cabo**
no modelo/contexto:

```python
context.external_radius_ground_return = np.array([...])   # (N_cabos,)
context.depth_ground_return           = np.array([...])   # (N_cabos,)
```

e `_cable_distance_matrices` + `earth_return_parameters` leem esses vetores.
Testável isoladamente; documenta explicitamente "qual raio o solo enxerga".

### 3.4 Efeito esperado

Config. 2, termo próprio de `Zg`/`Pg`: `x_jk` passa de `0.024` para `0.0578` m
→ `K0(γ_g · x)` menor, `ln`-termo maior → **`Zg` próprio menor, `Pg` próprio
maior** (capacitância de retorno menor). Combina com o `Yi` do duto (dielétrico
grande) para dar os `Z_cm`/`v_m` altos dos modos terra/entre-blindagens
(Figs 5.9–5.10) **sem** o `ε_a` do GMD.

### 3.5 Regressão

Config. 1 (`SingleCoreCableStrategy`, sem `enclosure`) — inalterada.
`hdpe_300mm2` / `hdpe_2000mm2` — passam a usar o raio do tubo; conferir contra
o MATLAB/COMSOL desses casos (pode **melhorar** a concordância).

---

## 4. Parte B — parâmetros internos a partir de FEM/COMSOL

### 4.1 Contrato de saída (igual ao analítico)

`InternalPerUnitParameters.matrices()` devolve hoje:

```python
{
  'impedance_matrix'            : Zi,   # (Nf, ΣM, ΣM)  complexo  [Ω/m]
  'shunt_admittance_matrix'     : Ye,   # (Nf, ΣM, ΣM)  complexo  [S/m]
  'potential_coefficient_matrix': Pi,   # (ΣM, ΣM)      real, freq-indep  [m/F]
  'capacitance_matrix'          : Ci,   # (ΣM, ΣM)      real, freq-indep  [F/m]
  'block_sizes'                 : [M]*N # condutores por cabo
}
```

De tudo isso, `quasi_tem_approx_matrices` só usa `impedance_matrix`,
`potential_coefficient_matrix` e `block_sizes`. O builder FEM precisa produzir
pelo menos esses três.

### 4.2 Novo builder

`analytical_forms/single_core_cable.py` (ou `utils/comsol_data.py`):

```python
class InternalParametersFromFEM:
    """Monta o dict `internal_matrices` (mesma forma que
    InternalPerUnitParameters.matrices()) a partir de matrizes por-cabo
    medidas em FEM/COMSOL. Para N cabos idênticos, tila kron(I_N, bloco).
    Para sistemas heterogêneos (SCC + ECC), monta bloco-diagonal, como
    _matrices_heterogeneous."""

    def __init__(self, frequencies, fem_cable_blocks):
        # fem_cable_blocks: lista de dicts, um por cabo físico:
        #   {'Zi': (Nf_fem, M, M) complexo [Ω/m],  'C': (M, M) real [F/m] nodal}
        # Zi é interpolado (log-log, Re/Im) na grade `frequencies`; C é constante.
        ...

    def matrices(self):
        # Zi_full  = block-diag / kron dos Zi(interp) por cabo
        # Pi_full  = block-diag / kron dos C^-1 por cabo
        # Ci_full  = block-diag / kron dos C
        # Ye_full  = jw * inv(Pi_full)         (só p/ inspeção)
        # block_sizes conforme M de cada cabo
        return {...}   # mesmas 6 chaves
```

### 4.3 Composição (inalterada)

```python
internal = InternalParametersFromFEM(freqs, blocks).matrices()
pul      = PerUnitParameters(mtl_hdpe, freqs)           # mtl com Parte A aplicada
earth    = pul.earth_return_parameters('magalhaes_xue', 'magalhaes_xue')
qt       = pul.quasi_tem_approx_matrices(internal, earth)
qt       = pul.propagation_matrices(qt)                 # γ_v, Zc, ...
```

---

## 5. Dados COMSOL — os arquivos existentes já servem

**Não é preciso nova exportação** e não há dependência de COMSOL neste plano.
Os dois arquivos em `testData/andreata_case2/Results/` já modelam a **geometria
excêntrica exata da Config. 2** (cabo em repouso no fundo do tubo HDPE, ar +
tubo na malha, fronteira externa na superfície do tubo — premissa confirmada).
`Zi` e `Yi` deles entram **integralmente** no caminho FEM-híbrido; o pyLCP já
tem os leitores. Resta só um método de leitura novo (§5.3) + interpolação de grade.

### 5.1 `cmsl_internal_impedance_matrix.txt` → `Zi(f)` 2×2

Magnetodinâmica FD, excitação por corrente. Leitor existente:
`ComsolPostProcessor.get_scc_internal_impedance_matrix_combined()` →
`{'scenarios': {'measured': {'impedance_matrix': Zi}}}` com `Zi` `(Nf, 2, 2)`
`[Ω/m]`, índices núcleo=0 / blindagem=1:

```
Zi[:,0,0] = data1(V_core_i)     núcleo excitado, tensão no núcleo   (próprio)
Zi[:,0,1] = Zi[:,1,0] = data1(V_sheath)   mútuo núcleo–blindagem
Zi[:,1,1] = data2(V_sheath_i)   blindagem excitada, tensão na blindagem (próprio)
```

Esse `Zi` **já contém `Ze`** — a indutância externa do anel ar + tubo até a
superfície externa do tubo (fronteira do modelo). É exatamente `Zi + Ze` da
decomposição de Andreata; a composição `Zs = Zi_fem + Zg` fecha com o retorno
pela terra analítico da Parte A (raio = `D2/2`).

### 5.2 `cmsl_internal_admittance_charge_method.txt` → `C` nodal 2×2

Eletrostática, método direto de carga, `εr` reais (ar = 1, HDPE = 2.35, XLPE
com correção semicondutora, PVC). Leitor existente:
`ComsolPostProcessor.get_scc_internal_admittance_matrix_combined()` → `Yi = jω·C`
com `C` **nodal** `(2, 2)` `[F/m]`, freq-independente:

```
C[0,0] = C_core_coreExc                    = Cc         (cap. núcleo–blindagem)
C[0,1] = C[1,0] = C_shIn_coreExc           = -Cc
C[1,1] = C_shIn_shExc + C_shOut_shExc      = Cc + Cs    (Cs = blindagem → sup. do tubo,
                                                          através de PVC + ar + HDPE)
```

Verificado numericamente: `Cs = C_shOut_shExc ≈ 1,37·10⁻¹⁰ F/m` — cerca de 12×
menor que a capacitância só do PVC (`~1,6·10⁻⁹`), coerente com o grande
dielétrico ar + HDPE em série. **O duto está capturado.**

Para o pipeline: `Pi_cabo = C⁻¹` (2×2), e `internal_matrices['potential_
coefficient_matrix'] = kron(I₃, Pi_cabo)` — mesma forma que o caminho analítico
(onde `capacitance_matrix` também é a matriz nodal, ver `compare_capacitance`).

### 5.3 Ações no pyLCP (só código, sem COMSOL)

- Novo `ComsolPostProcessor.get_scc_internal_capacitance_matrix_combined()` —
  devolve a matriz nodal `C` `(2, 2)` (ou `(7, 7)` p/ ECC) crua, sem o `jω`
  (hoje só existe a versão `·jω` para overlay).
- `InternalParametersFromFEM` consome `Zi` (5.1) + `C` (5.2), **interpola `Zi`
  na grade analítica** (log-log, Re e Im separados; `C` é constante) e tila
  `kron(I₃, ·)` → dict de 6 chaves.
- **Fronteira do modelo — confirmada:** a fronteira externa das duas `.mph` é a
  superfície externa do tubo HDPE (sem caixa de ar). `Zi_fem` = "tudo até
  `D2/2`" e casa direto com o `Zg` analítico da Parte A. Nenhuma verificação
  pendente.

### 5.4 Checagem cruzada (opcional)

`Le = Im{Zi_fem}/ω` em baixa frequência deve bater com `μ0 ε0 (C0)⁻¹` se um dia
houver um export `C0` (todos `εr = 1`) — Andreata eq. 4.48. Não bloqueia nada.

---

## 6. API do pipeline híbrido

`analytical_forms/single_core_cable.py` (ou novo `analytical_forms/hybrid.py`):

```python
def build_pul_matrices(mtl, frequencies, *,
                       internal_source='analytical',   # | 'fem'
                       fem_blocks=None,
                       zg_form='magalhaes_xue', yg_form='magalhaes_xue',
                       internal_form='hybrid'):
    """Orquestra interno (analítico OU FEM) + retorno pelo solo (analítico)
    + composição quasi-TEM + propagação de fase. Devolve o mesmo dict que
    quasi_tem_approx_matrices + propagation_matrices."""
    if internal_source == 'analytical':
        internal = InternalPerUnitParameters(mtl, frequencies).matrices(internal_form)
    elif internal_source == 'fem':
        internal = InternalParametersFromFEM(frequencies, fem_blocks).matrices()
    pul   = PerUnitParameters(mtl, frequencies)
    earth = pul.earth_return_parameters(zg_form, yg_form)
    qt    = pul.quasi_tem_approx_matrices(internal, earth)
    return pul.propagation_matrices(qt)
```

---

## 7. Integração — `andreata_case2` + modal

1. Novo cenário `'fem'` ao lado de `'1'/'2'/'3'`:

```python
fem_blocks = ComsolPostProcessor(__file__).get_hybrid_internal_blocks()   # 3× bloco 2×2
pul_data['scenarios']['fem'] = {
    'mtl': mtl_0,                      # flat_hdpe_enclosed_model (tubo real, Parte A ativa)
    'internal_source': 'fem', 'fem_blocks': fem_blocks,
    'zg_form': 'magalhaes_xue', 'yg_form': 'magalhaes_xue',
}
```

2. Loop de cenários: se `internal_source == 'fem'`, usar `build_pul_matrices`;
   senão, o caminho atual (GMD/ERS/underground fica para comparação).

3. **Decomposição modal roda sobre `'fem'`** (não mais `'3'`). Plotter
   `ModalPropagationPlotter` gera Figs 5.8/5.9/5.10 do cenário FEM;
   opcional: overlay tracejado do cenário `'3'` (GMD) para comparar.

4. `check_pul_passivity` no cenário `'fem'`.

---

## 8. Validação

| Alvo | Referência | Critério |
|---|---|---|
| `Zg`/`Pg` próprios com raio do tubo | cálculo manual `K0(γ_g·D2/2)` | exato |
| `Z'`, `Y'` FEM-híbrido | `andreata_case2_series_impedance_matrix.mat` / `..._shunt_admittance_matrix.mat` (referência do próprio Andreata, que **é** FEM) | erro rel. < ~5 % (molde: `validate_against_case1_reference`) |
| `Z'`, `Y'` FEM-híbrido vs GMD | cenário `'3'` | quantificar o quanto o GMD desvia |
| `α_m`, `v_m`, `\|Z_cm\|` FEM-híbrido | Figs 5.8–5.10 | forma + magnitude; **deve melhorar** vs GMD |
| passividade `Z'`, `Y'`, `Yc` FEM-híbrido | Gustavsen eq. 3 | PSD em toda a faixa |

---

## 9. Generalização (Configs 3–5)

Mesmo padrão:

| Config. | interno | raio externo p/ o solo |
|---|---|---|
| 3 (SCC + ECC enterrados) | analítico heterogêneo (já existe) | `max(r4_SCC, r_ECC)` por cabo — hoje OK (sem `enclosure`) |
| 4 (SCC em tubos + ECC compartilhando um tubo) | **FEM** (bloco 3×3 do tubo compartilhado + 2×2 dos outros) | tubo `D2/2` (Parte A já cobre `SingleCoreCableWithECCInHDPEStrategy`) |
| 5 (idem, ECC de cobre nu) | **FEM** (só muda o ECC na malha) | tubo `D2/2` |

`InternalParametersFromFEM` já previsto para blocos heterogêneos (bloco-diagonal
com `M` variável) → serve as três.

---

## 10. Fases

### Fase A — raio externo (0,5–1 dia)
- [ ] campo `external_radius_ground_return` / `depth_ground_return` nos strategies HDPE.
- [ ] `_cable_distance_matrices` (HDPE + shared-HDPE) lê o campo; termo próprio usa raio + profundidade do tubo.
- [ ] regressão: Config. 1 inalterada; `hdpe_300mm2`/`hdpe_2000mm2` conferidos.
- [ ] teste unitário: `diag(d_matrix_ground_return)` do `flat_hdpe_enclosed_model` = `[0.0578]*3`.

### Fase B — leitor + builder FEM (1 dia)
- [ ] `ComsolPostProcessor.get_scc_internal_capacitance_matrix_combined()` — matriz nodal `C` crua (sem `jω`).
- [ ] `ComsolPostProcessor.get_hybrid_internal_blocks()` — junta `Zi` (5.1) + `C` (5.2), devolve `[{'Zi','C'}]` por cabo.
- [ ] `InternalParametersFromFEM` — interpola `Zi` na grade analítica + tila `kron(I₃, ·)` → dict de 6 chaves. Prever bloco heterogêneo (ECC) desde já.

### Fase C — orquestrador + caso (1 dia)
- [ ] `build_pul_matrices(internal_source='analytical'|'fem', ...)`.
- [ ] `andreata_case2`: cenário `'fem'`, decomposição modal sobre ele, Figs 5.8–5.10 (+ overlay tracejado do GMD p/ comparação).
- [ ] `check_pul_passivity` + diagnósticos no cenário `'fem'`.

### Fase D — validação e doc (0,5–1 dia)
- [ ] tabela §8 preenchida (`Z'`/`Y'` vs `.mat` de Andreata; GMD vs FEM).
- [ ] atualizar `PLANO_CAP5_PROPAGACAO.md`, `DESENVOLVIMENTO_MODAL_CAP5.md`, README, `andreata_case2/README.md`.

**Total: ~3–4 dias.** Sem dependência de COMSOL — os dois arquivos já existem.

---

## 11. Riscos

1. ~~Fronteira externa da `.mph`~~ — **confirmado**: é a superfície externa do
   tubo, sem caixa de ar. `Zi_fem` e `Zg` casam em `D2/2` (§1, §5.3). Sem risco.
2. **Interface `Zi_fem` ↔ `Zg` no código.** A composição `Zs = Zi_fem + Zg` só é
   exata porque a fronteira do FEM e o `x_jk` do `Zg` (Parte A) são o mesmo
   círculo `D2/2`. Manter esse contrato documentado num único lugar
   (`InternalParametersFromFEM` + `_cable_distance_matrices`).
3. **Profundidade do objeto representativo.** Tubo centrado em `-1.1732`, cabo em
   `-1.2`. Definir e fixar (usar o centro do tubo quando o tubo representa o cabo).
4. **Grade de frequência do COMSOL ≠ grade analítica.** `get_general_parameters`
   já lida com isso para overlays; para *compor* é preciso interpolar o `Zi(f)`
   FEM na grade analítica (log-log, parte real e imaginária separadas).
5. **Camadas semicondutoras no FEM.** Ou modelar a geometria real, ou dobrar a
   permissividade como no lado analítico — escolher e casar com o `Zi`.
