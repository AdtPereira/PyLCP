# Desenvolvimento — Características de Propagação no Domínio Modal (Cap. 5 de Andreata)

**Data:** 2026-08.
**Escopo:** Configuração 1 da dissertação de Andreata — três cabos monopolares
(núcleo + blindagem) **diretamente enterrados** em arranjo plano, 34,5 kV.
6 condutores ⇒ 6 modos. Objetivo: a partir de `Z'(f)` e `Y'(f)` (já produzidos
pelo `andreata_case1`), calcular e plotar **`α_m`, `v_m`, `|Z_cm|`** por modo,
reproduzindo as **Figuras 5.5, 5.6 e 5.7**.

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
| `testData/andreata_case1/andreata_case1.py` | Chama `propagation_matrices`, roda `check_pul_passivity` + relatório, `ModalDecomposition` no cenário `p100_er1`, e `ModalPropagationPlotter.plot_all()` (3 figuras em `Results/`). |

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

### 3.8 Diagnóstico de passividade (`utils/passivity_check.py`)

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

## 5. Validação (`andreata_case1`, cenário `p100_er1`, solo ρ = 100 Ω·m constante)

```
mode labels                : (ground, inter_sheath_1, inter_sheath_2, coaxial_1, coaxial_2, coaxial_3)
max off-diag ratio Zm      : 2.3e-10          (≈ 0  → Z_m diagonal)
max scalar-relation error  : 8.4e-3           (Z_m,j·Y_m,j vs λ_j)
classification similarity  : 0.992  (min sobre os 6 modos)
reconstrução T_I Λ T_I⁻¹ ≈ Y'Z'               < 1e-16
passividade (Z', Y', Yc)   : PASSIVE  em toda a faixa
```

**Comparação com o PDF** (valores em 10⁶ Hz):

| grandeza | modo terra | entre-blindagens | coaxiais | Andreata (Figs 5.5–5.7) |
|---|---|---|---|---|
| `α_m` [Np/m] | 0,13 | 0,02–0,04 | 1,3·10⁻⁴ | terra ≫ blindagens ≫ coaxiais; cruzamento ~10² Hz ✓ |
| `v_m` [m/s] | 3,1·10⁷ | 3,5–4,1·10⁷ | 1,8·10⁸ | coaxiais ~2·10⁸, blindagens ~0,3–0,45·10⁸, terra ~0,2·10⁸ ✓ |
| `\|Z_cm\|` [Ω] | ~50–60 | ~16–20 | ~13 | terra ~40–60, blindagens ~15–25, coaxiais ~13 ✓ |

As 3 figuras batem com a dissertação em forma e magnitude, com desvio de
~10–20 % em alta frequência nos modos terra/entre-blindagens — **esperado**,
porque Andreata usa solo com parâmetros dependentes da frequência (Salvador
et al. 2020) e aqui o solo é constante. Fechar essa diferença é a Fase 3 do
plano (`models/frequency_dependent_soil.py`).

---

## 6. Limitações conhecidas

1. **Solo constante** (ρ = 100 Ω·m) em vez do modelo FD de Salvador (2020) —
   desvio de ~10–20 % em AF. Fase 3.
2. **Oscilação cosmética** dos 3 modos coaxiais (quase coincidentes) na última
   meia-década (5·10⁵–10⁷ Hz) por degenerescência genuína — o refino de cluster
   (§3.6) mitiga mas não elimina; próximo passo é ancorar o rastreamento nas
   **matrizes-limite assintóticas reais** de Wedepohl (1996) §7
   (BF: `eig(Re{C·R_dc})`; AF: `eig(Re{P⁻¹ Z'_ce})`).
3. **Só Configuração 1.** Configs 2–5 (duto HDPE, ECC, 7º modo, parâmetros
   externos por FEM) estão no plano, não implementadas.
4. `ModalPropagationPlotter` é autônomo (não usa `PLOT_CONFIG`); overlay de
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
