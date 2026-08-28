# Desenvolvimento — Características de Propagação no Domínio Modal (Cap. 5 de Andreata)

**Data:** 2026-08.
**Escopo:** Configurações **1 e 2** da dissertação de Andreata — três cabos
monopolares (núcleo + blindagem), arranjo plano, 34,5 kV; **diretamente
enterrados** (Config. 1) ou **em dutos HDPE individuais** (Config. 2).
6 condutores ⇒ 6 modos em ambas. Objetivo: a partir de `Z'(f)` e `Y'(f)` (já
produzidos por `andreata_case1` / `andreata_case2`), calcular e plotar
**`α_m`, `v_m`, `|Z_cm|`** por modo, reproduzindo as **Figuras 5.5–5.7**
(Config. 1) e **5.8–5.10** (Config. 2).

O plano completo (incl. fases futuras: Configs 2–5, solo dependente da
frequência) está em [`PLANO_CAP5_PROPAGACAO.md`](PLANO_CAP5_PROPAGACAO.md).

---

## 1. Arquivos

### Criados

| Arquivo | Conteúdo |
|---|---|
| `mtl_main/propagation.py` | `phase_domain_propagation(Zs, Ysh)` → `γ_v`, `γ_i`, `Z_c`, `Y_c` no domínio de fase (via `scipy.linalg.sqrtm`). Fatorado do laço que existia só em `overhead_lines.py`. |
| `analytical_forms/modal_analysis.py` | `ModalDecomposition`: `decompose()` → `track()` → `modal_parameters()` → `_classify()`. Helpers `default_scc_roles`, `canonical_permutation`, `REFERENCE_PATTERNS_6C`. |
| `plotter/modal_plotter.py` | `ModalPropagationPlotter` — Figs 5.5 / 5.6 / 5.7. `MODE_STYLE` (cor/traço fixos por rótulo de modo). |
| `utils/passivity_check.py` | `check_pul_passivity`, `print_passivity_report`, `min_real_part_eigenvalue` — **avaliação** de passividade (não enforcement). |

### Alterados

| Arquivo | Mudança |
|---|---|
| `analytical_forms/overhead_lines.py` | `PerUnitParameters.pul_matrices` passou a chamar `phase_domain_propagation` (resultado idêntico; regressão conferida em `ohtl_deConti_ex51`). |
| `analytical_forms/single_core_cable.py` | Novo `PerUnitParameters.propagation_matrices(quasi_tem_matrices)` — anexa `propagation_voltage_matrix`, `propagation_current_matrix`, `characteristic_impedance_matrix`, `characteristic_admittance_matrix` ao dict de `quasi_tem_approx_matrices` (retrocompatível). |
| `testData/andreata_case1/andreata_case1.py` | Chama `propagation_matrices`, roda `check_pul_passivity` + relatório, `ModalDecomposition` no cenário `p100_er1`, e `ModalPropagationPlotter.plot_all()` — Figs 5.5 / 5.6 / 5.7. |
| `testData/andreata_case2/andreata_case2.py` | Idem no cenário **`fem`** (pipeline FEM-híbrido — ver [`PLANO_PIPELINE_HIBRIDO.md`](PLANO_PIPELINE_HIBRIDO.md); `Zi`/`Yi` do COMSOL + retorno pela terra analítico com o raio do tubo, **sem** o GMD de Lafaia). Figs 5.8 / 5.9 / 5.10 + 4 gráficos de comparação vs. MATLAB (`self_impedance_phase_a_sheath`, `self_admittance_phase_a_sheath`, `earth_return_impedance_phase_a`, `earth_return_potential_coeff_phase_a`). |

---

## 2. Fundamentação (equações do Cap. 5, N = 6)

```
entrada:  Z'(f), Y'(f)  ∈ ℂ^{Nf×6×6}     (ordem canônica [c1 c2 c3 s1 s2 s3])

1. autodecomposição   λ, T_I  =  eig(Y' Z')                              (eq. 5.19)
2. matriz de tensão   T_V  =  (T_Iᵀ)⁻¹   ⟹  T_V⁻¹ = T_Iᵀ  (exato)         (eq. 5.20)
                      ⟹  Z' Y' = T_V · diag(λ) · T_V⁻¹                    (eq. 5.18)
3. rastreamento       reordena/refaseia colunas de T_I(f_k) para casar com
                      T_I(f_{k-1})  (switching-back — ver §3.5)
4. imped. série modal Z_m,j  =  [T_Iᵀ Z' T_I]_jj                          (eq. 5.21, diagonal)
5. const. de propag.  Γ_m  =  √λ ,  ramo físico  α_m = |Re|, β_m = |Im|   (eqs. 5.24, 5.30–5.31)
6. veloc. de fase     v_m  =  ω / β_m                                     (eq. 5.32)
7. imped. caract.     Z_cm  =  Z_m / Γ_m     Y_cm = 1 / Z_cm             (eqs. 5.23, 5.27)
8. rótulo do modo     projeta colunas de T_I(f_min) nos padrões           (§5.3, eqs. 5.34–5.36)
```

**Os 6 modos da Configuração 1** (rótulos usados no código):

| # | rótulo | descrição (Andreata §5.3) | vs f |
|---|---|---|---|
| 1 | `ground` | 3 blindagens em fase, retorno pela terra | ~invariante |
| 2 | `inter_sheath_1` | 1 blindagem × outra blindagem | ~invariante |
| 3 | `inter_sheath_2` | 1 blindagem × as outras duas | ~invariante |
| 4 | `coaxial_1` | 3 núcleos × terra (BF) → núcleo×blindagem (AF) | varia |
| 5 | `coaxial_2` | núcleo × núcleo (BF) → núcleo×blindagem (AF) | varia |
| 6 | `coaxial_3` | núcleo × 2 núcleos (BF) → núcleo×blindagem (AF) | varia |

---

## 3. Decisões de projeto (as não-óbvias)

### 3.1 Fatoração da propagação de fase (`mtl_main/propagation.py`)

O laço `sqrtm(Z'Y')` / `sqrtm(Y'Z')` / `Z_c = Y'⁻¹γ_i` / `Y_c = Z'⁻¹γ_v` existia
**apenas** em `analytical_forms/overhead_lines.py`; o caminho SCC parava em
`Z'/Y'`. Extraído para um módulo único, reusado por OHTL e SCC. `overhead_lines.py`
passou a chamá-lo sem alteração de resultado (mesma sequência LU + `sqrtm`).

### 3.2 `T_V = (T_Iᵀ)⁻¹` em vez de um segundo `eig`

`Z'Y' = (Y'Z')ᵀ` (pois `Z'`, `Y'` simétricas). Logo, se `Y'Z' = T_I Λ T_I⁻¹`,
então `Z'Y' = (T_I⁻¹)ᵀ Λ T_Iᵀ`, i.e. `T_V = (T_Iᵀ)⁻¹`. Isso:
- torna a eq. 5.20 (`T_V⁻¹ = T_Iᵀ`) **exata por construção** — sem tolerância;
- evita um segundo `eig(Z'Y')` e o problema de **casar** os autovalores/autovetores
  das duas decomposições.

### 3.3 `Z_m` via `diag(T_Iᵀ Z' T_I)`, **não** `Y_m`

Teorema (Wedepohl): `Z' t_j` é paralelo à j-ésima coluna de `T_V`, portanto
`Z_m = T_V⁻¹ Z' T_I = T_Iᵀ Z' T_I` é **diagonal** em aritmética exata.
`Y_m = T_I⁻¹ Y' T_V` também seria, mas depende dos **autovetores à esquerda**
(`T_I⁻¹`) e é mal condicionada perto de degenerescência. Medido no
`andreata_case1`: off-diagonal de `Z_m` ≈ 1e-10 (≈ 0), off-diagonal de `Y_m`
chegava a ~38 %. Por isso `Z_cm = Z_m / Γ_m` (e `Y_cm = 1/Z_cm`, domínio modal
é escalar desacoplado) — `Y_m` é devolvido só para inspeção.

### 3.4 `Γ_m = √λ` com ramo físico `α = |Re|`, `β = |Im|`

Para uma onda que propaga e atenua, `Γ = α + jβ` com `α, β > 0` ⟹
`Im(λ) = 2αβ > 0` ⟹ `√λ` principal já dá `Re > 0` e `Im > 0`. Mas perto de DC
o ruído numérico pode dar `Im(λ) < 0`, e `√` do plano próximo ao eixo real
negativo é hipersensível ao sinal de `Im` (`√(−1−0j) = −j`), gerando `v_m < 0`.
Solução: `α_m = |Re √λ|`, `β_m = |Im √λ|` (α e β são magnitudes). Documentado no
código.

### 3.5 Rastreamento = *switching-back procedure*

`ModalDecomposition.track()` implementa a mesma sub-rotina que Andreata (§5.2)
reusou: **Gustavsen (2008) §IV-A**, que é a correlação por produto escalar de
**Wedepohl-Nguyen-Irwin (1996) §6** (a função `intercheig` da *Matrix Fitting
Toolbox*). A cada passo casa as colunas de `T_I(f_k)` às de `T_I(f_{k-1})`
maximizando `Σ_j |⟨t_j(f_k), t_j(f_{k-1})⟩|`.

Diferenças em relação ao `intercheig` original (melhorias):
- atribuição resolvida de forma **ótima** (`scipy.optimize.linear_sum_assignment`,
  Hungarian) em vez do máximo-guloso linha-a-linha — mais robusto quando vários
  modos trocam de posição na mesma frequência (caso coaxial em AF);
- **âncora no meio da banda**, varrendo para os dois lados, em vez de varrer a
  partir de um extremo — na Config. 1 os modos coaxiais são mais separados nas
  décadas centrais.

Considerou-se e **não se adotou** o método de Newton-Raphson de Wedepohl (1996):
ele produz `T_I(ω)` *suave* e ajustável por funções racionais, o que só importa
para modelos de linha **no domínio modal** (EMTP). O pyLCP + NLT (Cap. 6 de
Andreata) é efetivamente **domínio de fase**; para plotar `α_m`, `v_m`, `|Z_cm|`
basta ordenação consistente dos modos, que o rastreamento por correlação
entrega. A própria Seção 6 do paper de Wedepohl usa esse método de correlação
como validação do NR ("agree perfectly").

### 3.6 Refino de clusters de autovalores degenerados

`decompose(cluster_rtol=3e-2)`: dentro de um cluster de autovalores que
concordam nessa tolerância relativa, `np.linalg.eig` devolve uma base
**arbitrária** do autoespaço. Reescolhe-se a base diagonalizando a impedância
série modal restrita `B = Uᵀ Z' U` (matriz complexa simétrica ⟹ seus
autovetores satisfazem `vᵢᵀ vⱼ = 0`, a mesma propriedade da eq. 7 de Gustavsen
2008). Sem isso, os 3 modos coaxiais quase-idênticos da Config. 1 geravam
oscilação em `|Z_cm|` na última meia-década (5·10⁵–10⁷ Hz).

### 3.7 Ordem canônica de condutores

O pyLCP monta `Z_i = kron(I₃, Z_ij)` ⟹ linhas `[c1 s1 c2 s2 c3 s3]`; Andreata
usa `[c1 c2 c3 s1 s2 s3]`. `canonical_permutation(default_scc_roles(...))`
converte (permutação `[0 2 4 1 3 5]` para o caso 3×SCC). Passividade e
autovalores são invariantes por essa permutação simétrica.

### 3.8 Classificação de modos em 2 estágios (`_classify`)

A primeira versão projetava as colunas de `T_I` **na frequência mínima** nos
padrões de referência (eqs. 5.34–5.36). Funcionou para a Config. 1
(similaridade 0,99) mas **inverteu** os rótulos na Config. 2: com o duto, em
BF os modos "de blindagem" carregam corrente de núcleo apreciável e os modos
só se separam acima de ~100 Hz.

Nova rotulagem, robusta para as duas configs (similaridade = 1,000):

1. **separa** os 3 modos de blindagem dos 3 coaxiais pela **fração de energia
   de `T_I` nas linhas dos núcleos**, média numa banda 100 Hz – 100 kHz (onde
   toda config está desacoplada). Os modos terra/entre-blindagens são
   invariantes com a frequência e não têm participação de núcleo (Andreata
   §5.3); os 3 com menor fração ⟹ grupo blindagem, os outros 3 ⟹ grupo coaxial;
2. **sub-rotula** cada grupo de 3 pelo padrão de linha correspondente
   (`[1,1,1]` = `ground`/`coaxial_1`; `[1,0,-1]` = `inter_sheath_1`/`coaxial_2`;
   `[-1,2,-1]` = `inter_sheath_2`/`coaxial_3`), via `linear_sum_assignment`,
   numa frequência média da banda.

Diagnósticos extras: `classification_ref_freq_hz`, `classification_core_fraction`.

### 3.9 Configuração 2 — duto HDPE (pipeline FEM-híbrido)

A estrutura de 6 modos é idêntica (Andreata eqs. 5.37–5.39 ≡ 5.34–5.36).
A decomposição modal roda sobre o **cenário `fem`** do `andreata_case2`
(ver [`PLANO_PIPELINE_HIBRIDO.md`](PLANO_PIPELINE_HIBRIDO.md)):

- **interno** (`Zi`, `Yi`): do COMSOL, geometria excêntrica exata com ar + tubo
  HDPE (`InternalParametersFromFEM`). Sem o truque GMD de permissividade
  equivalente de Lafaia.
- **retorno pela terra** (`Zg`, `Yg`): analítico (`magalhaes_xue`), com o termo
  próprio usando o **raio externo do tubo** (`D2/2`), não o do SCC
  (`_cable_external_geometry` em `mtl_main/strategy.py`).

Efeito físico: a corrente de retorno das blindagens vê ar + tubo HDPE
(dielétrico grande, `εr` baixo) ⟹ `Z_cm` e `v_m` dos modos terra/entre-blindagens
**bem maiores** que na Config. 1 (Andreata §5.4.1, Figs 5.8–5.10). O caminho
FEM-híbrido reproduz `Z'/Y'` da referência FEM de Andreata (`.mat`) a **< 2 %**
(vs. ~10–20 % do GMD), e a oscilação de AF dos modos coaxiais some.

### 3.10 Diagnóstico de passividade (`utils/passivity_check.py`)

Linha passiva ⟹ `eig(Re{Z'})≥0` (R' PSD), `eig(Re{Y'})≥0` (G' PSD),
`eig(Re{Y_c})≥0` (Gustavsen 2008, eq. 3) em toda a faixa. Autovalor negativo
além de ruído numérico aponta **problema de formulação** (retorno pela terra
com condutância negativa, `εr` de camada semicondutora errado, sinal trocado
numa aproximação). O módulo **só reporta** — reparar seria o pipeline
FRP/FMP + Hamiltoniana de Gustavsen (2008), deliberadamente fora do escopo:
aqui `Z'/Y'` são amostras de forma fechada, não macromodelo racional ajustado.
Tolerância relativa à escala de `|Re{M}|` por frequência (não acusa `G' ≈ 0` de
dielétrico sem perdas).

---

## 4. Uso

```python
from analytical_forms.modal_analysis import ModalDecomposition, default_scc_roles
from plotter.modal_plotter import ModalPropagationPlotter
from utils.passivity_check import check_pul_passivity, print_passivity_report

# qt = PerUnitParameters(...).quasi_tem_approx_matrices(...)
# qt = PerUnitParameters(...).propagation_matrices(qt)     # opcional: γ_v, Z_c, ...

print_passivity_report(check_pul_passivity(freq, qt))

md = ModalDecomposition(
    freq, qt['series_impedance_matrix'], qt['shunt_admittance_matrix'],
    conductor_roles=default_scc_roles(mtl.num_sc_cables, mtl.num_conductors_per_scc))
res = md.modal_parameters()           # dict: alpha, beta, vphase, Zcm, Ycm, TI, TV,
                                      #       mode_labels, diagnostics, ...
ModalPropagationPlotter(__file__, res, config_name='Configuracao 1').plot_all()
```

`res['diagnostics']`: `max_offdiag_ratio_Zm`, `max_scalar_relation_error`,
`reciprocity`, `passivity`, `classification_similarity`.

---

## 5. Validação (solo ρ = 100 Ω·m constante)

Diagnósticos das duas configs:

```
                             Config. 1 (case1, p100_er1)   Config. 2 (case2, cenário 3 = GMD)
mode labels                  (ground, inter_sheath_1/2, coaxial_1/2/3)   (idem)
max off-diag ratio Zm        2.3e-10                       5.5e-12          (≈ 0 → Z_m diagonal)
max scalar-relation error    8.4e-3                        8.4e-3           (Z_m,j·Y_m,j vs λ_j)
classification similarity    1.000                         1.000            (min sobre os 6 modos)
reconstrução T_I Λ T_I⁻¹     < 1e-16                        < 1e-16
passividade (Z', Y', Yc)     PASSIVE                        PASSIVE
```

**Comparação com o PDF** (valores em 10⁶ Hz):

*Config. 1 — Figs 5.5–5.7:*

| grandeza | modo terra | entre-blindagens | coaxiais | Andreata |
|---|---|---|---|---|
| `α_m` [Np/m] | 0,13 | 0,02–0,04 | 1,3·10⁻⁴ | terra ≫ blindagens ≫ coaxiais; cruzamento ~10² Hz ✓ |
| `v_m` [m/s] | 3,1·10⁷ | 3,5–4,1·10⁷ | 1,8·10⁸ | ✓ |
| `\|Z_cm\|` [Ω] | ~50–60 | ~16–20 | ~13 | terra ~40–60, blindagens ~15–25, coaxiais ~13 ✓ |

*Config. 2 — Figs 5.8–5.10:*

| grandeza | modo terra | entre-blindagens | coaxiais | Andreata |
|---|---|---|---|---|
| `α_m` [Np/m] | ~0,1 | ~0,02–0,04 | ~2·10⁻⁴ | terra ≫ blindagens ≫ coaxiais ✓ |
| `v_m` [m/s] | ~5·10⁷ | **1,2–1,45·10⁸** | ~1,8·10⁸ | entre-blindagens **muito** mais rápidos que na Config. 1 ✓ |
| `\|Z_cm\|` [Ω] | ~130–200 | **~55–70** | ~13 | terra ~130–150, blindagens ~60–70, coaxiais ~13 ✓ |

A assinatura da Config. 2 (Andreata §5.4.1) — `Z_cm` e `v_m` dos modos
terra/entre-blindagens **bem maiores** que na Config. 1 por causa da isolação
ar + tubo HDPE — é reproduzida.

**Validação `Z'`/`Y'`/`Zg`/`Pg` vs. MATLAB** (referência FEM de Andreata, cenário
`fem` do `andreata_case2` — 4 gráficos de comparação em `Results/`):

| grandeza | erro rel. máx. |
|---|---|
| `Z'` série completa | 0,40 % |
| `Y'` shunt completa | 1,85 % |
| `Zg` retorno pela terra | 0,07 % |
| `Pg` coef. de potencial | 0,81 % |

---

## 6. Limitações conhecidas

1. **Solo constante** (ρ = 100 Ω·m) em vez do modelo FD de Salvador (2020) —
   afeta os modos terra/entre-blindagens em alta f. Fase 3 (a comparação
   `Z'/Y'` vs. MATLAB acima, < 2 %, sugere que o efeito é menor do que o
   estimado antes).
2. **Oscilação cosmética** dos 3 modos coaxiais na Config. 1 (grade até 10 MHz):
   degenerescência genuína — o refino de cluster (§3.6) mitiga. **Na Config. 2
   com o `Zi` FEM a oscilação some** (o `Zi` FEM tem a indutância externa real,
   suave). Próximo passo p/ a Config. 1: ancorar o rastreamento nas
   matrizes-limite assintóticas de Wedepohl (1996) §7.
3. **Offset de ~2 % em `C₂₂`** da Config. 2: diferença entre o FEM do COMSOL
   (usado pelo pyLCP) e o FEM próprio de Andreata — não é erro do pipeline.
4. **Configs 3–5 não implementadas** (ECC, 7º modo, duto compartilhado, família
   §5.4.2). No plano; a Parte A do pipeline híbrido já cobre a geometria de 4–5.
5. `ModalPropagationPlotter` é autônomo (não usa `PLOT_CONFIG`); overlay de
   referência MATLAB modal ainda não implementado.

---

## 7. Referências

- **Andreata, L. E. B.** (2025). *Análise das Características de Propagação e de
  Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não
  Metálicos no Contexto de Parques Eólicos.* Dissertação, PPGEE/UFMG. Cap. 5.
- **Wedepohl, L. M.; Nguyen, H. V.; Irwin, G. D.** (1996). Frequency-Dependent
  Transformation Matrices for Untransposed Transmission Lines using
  Newton-Raphson Method. *IEEE Trans. Power Delivery*, 11(3), 1538–1546.
  (§6 — correlação de autovetores; §7 — assintóticas.)
- **Gustavsen, B.** (2008). Fast Passivity Enforcement for Pole-Residue Models
  by Perturbation of Residue Matrix Eigenvalues. *IEEE Trans. Power Delivery*,
  23(4), 2278–2285. (§IV-A — switching-back; eq. 3 — critério de passividade;
  eq. 7 — perturbação de autovalor / ortogonalidade `vᵀv`.)
- **Paul, C. R.** (2008). *Analysis of Multiconductor Transmission Lines*, 2ª ed.
  (teoria modal — Cap. 7.)
