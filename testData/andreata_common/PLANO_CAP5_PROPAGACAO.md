# Plano — Características de Propagação no Domínio Modal (Cap. 5 de Andreata) — **somente Configuração 1**

---

## STATUS (atualizado)

**Fases 0–4 do núcleo implementadas e validadas contra as Figs 5.5 / 5.6 / 5.7.**

| Entregue | Arquivo |
|---|---|
| `phase_domain_propagation(Zs, Ysh)` (γ_v, γ_i, Zc, Yc) — fatorado do OHTL | `mtl_main/propagation.py` |
| `overhead_lines.py` migrado p/ o módulo compartilhado (regressão OK: `ohtl_deConti_ex51`) | `analytical_forms/overhead_lines.py` |
| `PerUnitParameters.propagation_matrices()` p/ SCC | `analytical_forms/single_core_cable.py` |
| `ModalDecomposition` — decompose (com refino de clusters degenerados) → track (âncora no meio da banda, varre p/ os dois lados) → modal_parameters → classify (6 modos, eqs 5.34–5.36) | `analytical_forms/modal_analysis.py` |
| `ModalPropagationPlotter` — Figs 5.5 (α_m), 5.6 (v_m), 5.7 (\|Z_cm\|) | `plotter/modal_plotter.py` |
| `andreata_case1.py` chama decomposição modal + gera as 3 figuras em `Results/` | `testData/andreata_case1/andreata_case1.py` |
| **Diagnóstico de passividade** (só avaliação — Gustavsen 2008 eq. 3): `eig(Re{Z'})`, `eig(Re{Y'})`, `eig(Re{Yc})` ≥ 0 ao longo da varredura | `utils/passivity_check.py` |
| Rastreamento de modos citado como *switching-back procedure* (Gustavsen 2008 §IV-A; Wedepohl 1996 §6 — sub-rotina `intercheig` que Andreata §5.2 reusou); atribuição ótima via Hungarian em vez do máximo-guloso de Gustavsen | docstrings de `modal_analysis.py` |

**Resultados** (cenário `p100_er1`, solo ρ=100 Ωm constante): rótulos limpos
`(ground, inter_sheath_1, inter_sheath_2, coaxial_1, coaxial_2, coaxial_3)`,
similaridade de classificação ≥ 0,99; `max_offdiag_ratio_Zm ≈ 0`,
`max_scalar_relation_error ≈ 1 %`; reconstrução `T_I Λ T_I⁻¹ ≈ Y'Z'` < 1e-16.
Passividade: `Z'`, `Y'`, `Yc` todos PSD em toda a faixa (`overall: PASSIVE`).
As 3 figuras batem com o PDF em forma e magnitude (dentro de ~10–20 %, coerente
com o uso de solo constante em vez do modelo FD de Salvador 2020).

**Pendências / próximos passos**
- Fase 3 (opcional): solo dependente da frequência (`models/frequency_dependent_soil.py`) — fecha o gap de ~10–20 % em alta f nos modos terra/entre-blindagens.
- Cosmético: pequena oscilação/queda dos 3 modos coaxiais (quase coincidentes) na última meia-década (5e5–1e7 Hz) por degenerescência — risco §6.1; mitigável com colapso de `Z_cm` ao valor médio dentro de cluster totalmente degenerado.
- `plot_config.py`: hoje o `ModalPropagationPlotter` é autônomo (não usa `PLOT_CONFIG`); adicionar entradas se quiser sobrepor referência MATLAB modal.
- Validação quantitativa figura-a-figura (tabela de erro por modo em f ∈ {1e2, 1e4, 1e6} Hz).
- `f_c` (eq. 5.33) como anotação nos gráficos.

---


> Escopo restrito: **Configuração 1** (Fig. 5.1) — três cabos monopolares (núcleo + blindagem)
> **diretamente enterrados** no solo, arranjo plano, 34,5 kV. São **6 condutores → 6 modos**.
> Objetivo: a partir de `Z'(f)` e `Y'(f)` (já produzidos pelo `andreata_case1`), calcular e plotar
> **`α_m`, `v_m`, `|Z_Cm|`** por modo, reproduzindo as **Figuras 5.5, 5.6 e 5.7** da dissertação.
>
> Fora de escopo (adiado): Configurações 2–5, condutor de aterramento (ECC) e seu 7º modo,
> duto HDPE / parâmetros externos por FEM, família de gráficos "por modo entre configurações"
> (§5.4.2, Figs 5.20–5.31).

---

## 1. Diagnóstico — o que já existe vs. o que falta (Config. 1)

### 1.1 Já implementado (reaproveitar)

| Item | Onde |
|---|---|
| `Z'(f)`, `Y'(f)` completos (Zi+Ze+Zg / Yi ∥ Yg) p/ 3 SCC enterrados em arranjo plano | `analytical_forms/single_core_cable.py::PerUnitParameters.quasi_tem_approx_matrices` → chaves `series_impedance_matrix`, `shunt_admittance_matrix` |
| Impedância/admitância interna (Bessel/aprox/híbrido) + correção de camada semicondutora | `InternalPerUnitParameters`, `apply_semiconducting_layer_correction` |
| Retorno pelo solo quasi-TEM (De Conti/Duarte/Alipio 2023 ≈ `magalhaes_xue`) | `PerUnitParameters.earth_return_parameters` |
| Caso montado até `Z'/Y'` com o solo `ρ = 100 Ωm` (cenário `p100_er1`) | `testData/andreata_case1/andreata_case1.py` |
| Referência MATLAB de `Z'/Y'` + remapeamento de ordem de condutor (`conductor_order=[0,3,1,4,2,5]`) | `utils/matlab_data.py::MatlabDataReader.get_scc_scenario_data` |
| **Propagação de FASE** `γ_v = sqrtm(Z'Y')`, `γ_i`, `Z_c`, `Y_c` (loop por frequência) | **somente** `analytical_forms/overhead_lines.py::PerUnitParameters.pul_matrices` (linhas ~516–558) |
| Plot atenuação `Re γ` × velocidade `ω/Im γ` p/ 1 elemento `[p,q]` | `plotter/ohtl_models.py::_propagation_subplots`, `plotter/lima_models.py` |
| Esquemático de seção transversal | `mtl_main/graphics.py::GroundReturnMTLRepresentation` |

### 1.2 Lacunas a preencher

1. **Transformação modal ausente.** Não há autovalores/autovetores de `Y'Z'` / `Z'Y'`, nem `T_I`, `T_V`.
2. **Propagação de fase não existe para SCC.** `single_core_cable.py::PerUnitParameters` para em `Z'/Y'`; só o caminho OHTL calcula `γ`, `Z_c`.
3. **Rastreamento de modos (mode tracking).** Autovalores trocam de posição entre frequências vizinhas (Andreata §5.2). Falta a continuidade dos autovetores (produto escalar frequência‑a‑frequência — Gustavsen 2008).
4. **Parâmetros modais.** `Γ_m = √λ` (ramo `Re ≥ 0`), `α_m = Re Γ_m`, `β_m = Im Γ_m`, `v_m = ω/Im Γ_m`, `Z_Cm = Γ_m/Y_m` (Andreata eqs. 5.23–5.32).
5. **Identificação/rotulagem dos 6 modos** (terra; entre‑blindagens 1 e 2; coaxiais 1/2/3) por inspeção das colunas de `T_I` (Andreata eqs. 5.34–5.36).
6. **Plots do domínio modal.** 3 figuras (`α_m`, `v_m`, `|Z_Cm|`) com os 6 modos sobrepostos → Figs 5.5, 5.6, 5.7.
7. **Solo com parâmetros dependentes da frequência** (Salvador et al. 2020, `ρ_LF = 100 Ωm`). Hoje `σ_g`, `ε_g` são escalares fixos. *(opcional — Fase 3; a versão inicial usa `ρ = 100 Ωm` constante como no `andreata_case1`)*
8. **Frequência de corte `f_c`** de desacoplamento (eq. 5.33) — utilitário informativo p/ anotar gráficos. *(opcional)*

---

## 2. Fundamentos a implementar (equações do Cap. 5, N = 6)

```
entrada:  Z'(f) [Nf,6,6] , Y'(f) [Nf,6,6]

1. produtos              A_i = Y'Z'   ,  A_v = Z'Y'                        (eqs. 5.14–5.15)
2. autodecomposição      A_i = T_I · diag(λ) · T_I^{-1}                    (eq. 5.19)
                         A_v = T_V · diag(λ) · T_V^{-1}   (mesmos λ)       (eq. 5.18)
                         validar  T_V^{-1} ≈ T_I^{t}                       (eq. 5.20)
3. rastreamento          reordenar/refasear colunas de T_I(f_k) p/ casar com T_I(f_{k-1})
                         (produtos escalares; semear da alta frequência)
4. constante modal       Γ_m(f) = √λ(f) ,  ramo  Re{Γ_m} ≥ 0             (eq. 5.24)
5. matrizes modais       Z_m = T_V^{-1} Z' T_I  (diagonal)                (eq. 5.21)
                         Y_m = T_I^{-1} Y' T_V  (diagonal)                (eq. 5.22)
6. imped. caract. modal  Z_Cm = Γ_m / Y_m                                 (eqs. 5.23, 5.27)
7. saída                 α_m = Re Γ_m   β_m = Im Γ_m   v_m = ω / Im Γ_m   (eqs. 5.30–5.32)
8. rótulo do modo        classify(T_I, papéis core/sheath)                (§5.3, eqs. 5.34–5.36)
```

**Os 6 modos da Configuração 1 (Andreata §5.3):**

| # | Rótulo | Descrição | Comportamento vs f |
|---|---|---|---|
| 1 | `ground` | corrente pelas 3 blindagens, retorno pela terra | ~invariante |
| 2 | `inter_sheath_1` | 1 blindagem × outra blindagem | ~invariante |
| 3 | `inter_sheath_2` | 1 blindagem × as outras duas | ~invariante |
| 4 | `coaxial_1` | 3 núcleos × terra (BF) → núcleo×blindagem (AF) | varia com f |
| 5 | `coaxial_2` | núcleo × núcleo (BF) → núcleo×blindagem (AF) | varia com f |
| 6 | `coaxial_3` | núcleo × 2 núcleos (BF) → núcleo×blindagem (AF) | varia com f |

`f_c = ρ_sh / [ π μ_sh (r_sh,out − r_sh,in)^2 ]`  (eq. 5.33) — opcional, anota o início do desacoplamento.

---

## 3. Arquitetura proposta

### 3.1 Novos artefatos (enxuto)

```
mtl_main/propagation.py              # phase-domain: γ_v, γ_i, Zc, Yc (extrai o loop de overhead_lines.py)
analytical_forms/modal_analysis.py   # ModalDecomposition: eig, tracking, params modais + classify (6 modos)
plotter/modal_models.py              # config declarativa dos 3 gráficos modais (estilo ohtl_models.py)
plotter/modal_plotter.py             # ModalPropagationPlotter(BasePlotter)
models/frequency_dependent_soil.py   # σ_g(f), ε_g(f) — Fase 3, opcional
```

Alterações:
```
analytical_forms/overhead_lines.py       # usar mtl_main/propagation.py (sem mudar resultado)
analytical_forms/single_core_cable.py    # PerUnitParameters.propagation_matrices();
                                         # (Fase 3) aceitar σ_g/ε_g vetoriais
testData/andreata_case1/andreata_case1.py    # chamar decomposição modal + novo plotter
testData/andreata_case1/plot_config.py       # 3 entradas: modal_attenuation / modal_phase_velocity / modal_char_impedance
```

### 3.2 `mtl_main/propagation.py` (domínio de fase)

```python
def phase_domain_propagation(Zs, Ysh):
    """Zs, Ysh: (Nf, N, N). Retorna:
       propagation_voltage_matrix   γ_v = sqrtm(Zs·Ysh)
       propagation_current_matrix   γ_i = sqrtm(Ysh·Zs)
       characteristic_impedance_matrix  Zc = Ysh^{-1}·γ_i
       characteristic_admittance_matrix Yc = Zs^{-1}·γ_v
    """
```
- Substituir o bloco duplicado em `overhead_lines.py::pul_matrices` por chamada a esta função;
  rodar `run_ohtl.py` e conferir que as figuras de Lima não mudam (regressão).
- Em `single_core_cable.py::PerUnitParameters`, acrescentar
  `propagation_matrices(quasi_tem_matrices)` → devolve as 4 chaves acima anexadas ao dict
  (retrocompatível).

### 3.3 `analytical_forms/modal_analysis.py`

```python
CANONICAL_ORDER = 'cores_then_sheaths'   # [c1,c2,c3,s1,s2,s3]

class ModalDecomposition:
    def __init__(self, freq, Zs, Ysh, conductor_roles):
        # conductor_roles: lista alinhada às linhas de Zs/Ysh, ex.
        #   [('core',0),('sheath',0),('core',1),('sheath',1),('core',2),('sheath',2)]
        # deriva a permutação P para a ordem canônica [cores | sheaths]
        ...

    def decompose(self):
        # p/ cada f (na ordem canônica):
        #   λ, TI = eig(Ysh @ Zs)
        #   _,  TV = eig(Zs @ Ysh)
        #   casar autovalores de TV aos de TI (nearest-λ) e reordenar TV
        #   normalizar autovetores: ∞-norm = 1; fase da 1ª componente não-nula = 0
        # retorna raw (sem tracking)

    def track(self):
        # continuidade por produto escalar; semente = frequência MAIS ALTA (modos ~ desacoplados)
        # varrer p/ baixo: no f_k, achar permutação+sinais que maximizam
        #   Σ_j |< TI_j(f_k) , conj(TI_j(f_{k+1})) >|   via scipy.optimize.linear_sum_assignment
        #   sobre a matriz de similaridade S = |TI(f_k)^H · TI(f_{k+1})|
        # aplicar a mesma permutação/sinais a λ, TI, TV

    def modal_parameters(self):
        Gamma = np.sqrt(lam);  Gamma *= np.where(Gamma.real < 0, -1, 1)   # Re ≥ 0
        Zm = einsum('fij,fjk,fkl->fil', inv(TV), Zs, TI).diagonal(axis1=1, axis2=2)
        Ym = einsum('fij,fjk,fkl->fil', inv(TI), Ysh, TV).diagonal(axis1=1, axis2=2)
        Zcm = Gamma / Ym
        return dict(gamma=Gamma, alpha=Gamma.real, beta=Gamma.imag,
                    vphase=self.w[:,None]/Gamma.imag, Zcm=Zcm, Ycm=Gamma/Zm,
                    TI=TI, TV=TV, lam=lam, mode_labels=self.classify(),
                    decoupling_freq=self._fc())   # arrays (Nf, 6)

    def classify(self, ref='high'):
        # projeta as colunas de TI (em f alta, e concilia com f média/baixa) sobre
        # REFERENCE_PATTERNS_6C (abaixo); rótulo = maior |projeção|
```

`REFERENCE_PATTERNS_6C` — colunas de `T_I` normalizadas na ordem canônica `[c1,c2,c3,s1,s2,s3]`,
transcritas das eqs. 5.34–5.36 de Andreata:

```python
REFERENCE_PATTERNS_6C = {
    # modos "de blindagem" (colunas 1-3 de T_I, ~invariantes com f)
    'ground'        : [ 0,   0,   0,   1/3, 1/3, 1/3],   # sheaths co-fase, retorno pela terra
    'inter_sheath_1': [ 0,   0,   0,   1/2,   0, -1/2],
    'inter_sheath_2': [ 0,   0,   0,  -1/3, 2/3, -1/3],
    # modos coaxiais (colunas 4-6, aqui no limite de ALTA frequência: núcleo × própria blindagem)
    'coaxial_1'     : [ 1,   0,   0,  -1,   0,   0 ],
    'coaxial_2'     : [ 0,   1,   0,   0,  -1,   0 ],
    'coaxial_3'     : [ 0,   0,   1,   0,   0,  -1 ],
}
```
> Ajustar os sinais/valores exatos conforme a convenção de montagem de `Z'/Y'` do pyLCP
> (validar contra as matrizes impressas por Andreata nas eqs. 5.34–5.36). A classificação
> usa **projeção normalizada** (robusta a escala/sinal), não igualdade.

Notas:
- **Reciprocidade**: conferir `‖T_V^{-1} − T_I^{T}‖ / ‖T_I‖` (aviso, não erro).
- **Degenerescência** dos modos coaxiais em BF: o `linear_sum_assignment` sobre a matriz de
  similaridade resolve a troca de rótulo; se persistir, casar no subespaço (projeção do bloco).
- Reaproveitar rotina de perturbação de autovalores de Gustavsen (2008) se já existir algo
  em `mom_so/` — checar antes de reimplementar.

### 3.4 Ordem canônica de condutores

O pyLCP monta `Z_i = kron(I_3, Z_ij)` → linhas `[c1,s1,c2,s2,c3,s3]`. Andreata usa
`[c1,c2,c3,s1,s2,s3]`. Definir **uma** matriz de permutação `P` (reutilizar a lógica de
`conductor_order=[0,3,1,4,2,5]` já usada p/ o MATLAB), aplicá-la a `Z'` e `Y'` **antes** da
decomposição, e documentá-la num único lugar (`modal_analysis.py`).

### 3.5 `models/frequency_dependent_soil.py` *(Fase 3 — opcional)*

```python
def soil_parameters(freq, rho_lf, model='alipio'):
    """Retorna sigma_g(f) [S/m], eps_r_g(f). Modelos: 'alipio' (Alipio-Visacro 2014),
       'salvador' (Salvador et al. 2020, usado por Andreata §5.4), 'constant'."""
```
Integração: `PerUnitParameters.__init__` já constrói `k_earth2`, `gamma_earth`, `jw_2pi_sg`
como vetores `(Nf,)` — basta permitir `self.sigma_1`, `self.e1` como arrays e ajustar broadcasting.
Validar isoladamente (modo terra e entre‑blindagens em AF são os mais sensíveis).

### 3.6 Plotters

`plotter/modal_models.py` — config declarativa (estilo `ohtl_models.py`), 3 gráficos, série = 1 modo:
```python
MODE_STYLE = {                       # cores/traços fixos por rótulo (consistentes nas 3 figs)
  'ground'        : {'color':'black','linestyle':'-'},
  'inter_sheath_1': {'color':'tab:blue','linestyle':'--'},
  'inter_sheath_2': {'color':'tab:blue','linestyle':'-.'},
  'coaxial_1'     : {'color':'tab:red','linestyle':'-'},
  'coaxial_2'     : {'color':'tab:green','linestyle':'-'},
  'coaxial_3'     : {'color':'tab:orange','linestyle':':'},
}
```

`plotter/modal_plotter.py::ModalPropagationPlotter(BasePlotter)`:
```python
def modal_attenuation(self):      # Fig. 5.5  — α_m (Np/m) × f, log-log, 6 curvas
def modal_phase_velocity(self):   # Fig. 5.6  — v_m (m/s, absoluto) × f, semilog-x, 6 curvas
def modal_char_impedance(self):   # Fig. 5.7  — |Z_Cm| (Ω) × f, semilog-x, 6 curvas
```
- unidades como Andreata: `α_m` em **Np/m** (parametrizável Np/km), `v_m` em **m/s absoluto**
  (não normalizado), `|Z_Cm|` em **Ω**; eixo x `1e-1 … 1e6 Hz`.
- reutilizar `save_figure`, `_format_axis`; overlay opcional de referência MATLAB modal (Fase 4).

---

## 4. Fases de implementação

### Fase 0 — Preparação (~0,5 dia)
- [ ] Fatorar `phase_domain_propagation` em `mtl_main/propagation.py`; migrar `overhead_lines.py`; rodar `run_ohtl.py` (regressão: figuras de Lima idênticas).
- [ ] Definir a matriz de permutação canônica `[cores | sheaths]` em `modal_analysis.py`.

### Fase 1 — Propagação de fase p/ SCC + validação (~1 dia)
- [ ] `PerUnitParameters.propagation_matrices()` em `single_core_cable.py`.
- [ ] Anexar `γ_v, γ_i, Zc, Yc` ao dict do cenário `p100_er1` no `andreata_case1.py`.
- [ ] Sanidade: `Z'`, `Y'` simétricas; `‖γ_v − γ_iᵀ‖` pequeno; `Re γ_v[0,0] ≥ 0`.

### Fase 2 — Decomposição e rastreamento modal (~1,5–2 dias)
- [ ] `ModalDecomposition.decompose()` + `modal_parameters()` (sem tracking).
- [ ] `track()` (produto escalar + `linear_sum_assignment`, semente na alta frequência).
- [ ] Testes:
  - reconstrução `T_I diag(λ) T_I^{-1} ≈ Y'Z'` (erro rel. < 1e-8);
  - `Γ_m² ≈ autovalores(Z'Y')` (< 1e-8);
  - `α_m ≥ 0` em toda a faixa; `v_m ≤ c` (com folga p/ ruído em BF);
  - colunas `ground` / `inter_sheath_*` de `T_I` ~constantes vs f.
- [ ] Reproduzir **Fig. 5.5 (α_m)** e **Fig. 5.6 (v_m)** — 6 curvas.

### Fase 3 — Identificação de modos + `|Z_Cm|` (~1 dia)
- [ ] `classify()` com `REFERENCE_PATTERNS_6C` (eqs. 5.34–5.36); conferir rótulos contra o texto de Andreata §5.3 (modos 1–3 ~invariantes; 4–6 coaxiais).
- [ ] Reproduzir **Fig. 5.7 (|Z_Cm|)**.
- [ ] *(opcional)* `frequency_dependent_soil.py` (Salvador 2020, `ρ_LF = 100 Ωm`) e re‑rodar as 3 figuras.
- [ ] *(opcional)* anotar `f_c` (eq. 5.33) nos gráficos.

### Fase 4 — Plots, driver e validação (~1 dia)
- [ ] `modal_models.py` + `ModalPropagationPlotter`; 3 entradas em `andreata_case1/plot_config.py`.
- [ ] `andreata_case1.py` chama os 3 plots; salva em `testData/andreata_case1/Results/`.
- [ ] Comparação figura‑a‑figura com o PDF: tabela de erro relativo por modo em `f ∈ {1e2, 1e4, 1e6} Hz`.
- [ ] *(opcional)* ingestão de referência MATLAB modal (`alpha_m`, `vphase_m`, `Zcm`) como scatter.
- [ ] Atualizar `README.md` (seção Cap. 5, Config. 1) e `testData/andreata_case1/BUGS_AND_FIXES.md`.

**Estimativa total: ~4–6 dias.**

---

## 5. Validação

| Nível | Referência | Critério |
|---|---|---|
| Unitário | reconstrução `T_I Λ T_I^{-1}` vs `Y'Z'` | erro rel. < 1e-8 |
| Unitário | `Γ_m²` vs autovalores de `Z'Y'` | erro rel. < 1e-8 |
| Físico | `α_m ≥ 0`; `v_m → c/√ε_eff` em alta f p/ modos coaxiais | qualitativo |
| Cruzado | `Z'`, `Y'` (cenário `p100_er1`) vs dump MATLAB do `andreata_case1` | reusar molde de `validate_against_case3_reference` |
| Cruzado | modo `coaxial_*` em AF vs `γ` de fase do núcleo isolado (`scc_34kV_andreata`) | < 5 % |
| Figuras | Figs 5.5 / 5.6 / 5.7 do PDF | inspeção visual + erro rel. por modo em 3 décadas de f |
| Rótulos | descrição textual §5.3 de Andreata | correspondência 1:1 dos 6 nomes de modo |

---

## 6. Riscos e limitações

1. **Mode tracking em degenerescência** (modos coaxiais quase idênticos em BF) pode trocar rótulo;
   mitigar com semeadura em alta frequência + `linear_sum_assignment` na matriz de similaridade.
2. **Ordem de condutores**: inconsistência entre `[c1,s1,…]` (pyLCP) e `[c1,c2,c3,s1,…]` (Andreata)
   quebra a identificação — centralizar numa única permutação e testá-la.
3. **Solo constante vs FD**: a versão inicial (`ρ = 100 Ωm` constante) desvia das Figs 5.5–5.7 em
   alta f (modo terra e entre‑blindagens); Fase 3 fecha a lacuna.
4. **`sqrtm` vs autovalores**: usar a via de autovalores como fonte de verdade p/ `Γ_m`; `sqrtm`
   só p/ conferência do caminho de fase.
5. **Camadas semicondutoras**: já tratadas via `ε_ins_eq` em `apply_semiconducting_layer_correction`;
   confirmar que a correção é aplicada **antes** da decomposição modal.

---

## 7. Checklist de arquivos

**Criar**
- `mtl_main/propagation.py`
- `analytical_forms/modal_analysis.py`
- `plotter/modal_models.py`, `plotter/modal_plotter.py`
- `models/frequency_dependent_soil.py`  *(Fase 3, opcional)*

**Alterar**
- `analytical_forms/overhead_lines.py` — usar `propagation.py` (sem mudar resultado)
- `analytical_forms/single_core_cable.py` — `PerUnitParameters.propagation_matrices()`; *(Fase 3)* `σ_g(f)`, `ε_g(f)` vetoriais
- `testData/andreata_case1/andreata_case1.py` — decomposição modal + 3 plots
- `testData/andreata_case1/plot_config.py` — `modal_attenuation` / `modal_phase_velocity` / `modal_char_impedance`
- `README.md`, `testData/andreata_case1/BUGS_AND_FIXES.md`
