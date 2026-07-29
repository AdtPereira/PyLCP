# Diagnósticos e correções — `andreata_case3` (Configuração 3, Figura 5.3)

**Contexto:** `andreata_case3` modela a Configuração 3 da referência (Figura 5.3): três
cabos SCC de potência (núcleo + blindagem) em arranjo plano, diretamente enterrados no
solo — **sem duto HDPE** —, mais um cabo de aterramento isolado (ECC) próximo ao cabo
mais à direita, sem encostar nele. É o primeiro caso do repositório com uma mistura
heterogênea de cabos (3 SCC de 2 condutores cada + 1 ECC de 1 condutor), o que expôs
uma série de suposições de "N cabos idênticos" implícitas em várias partes do pipeline.

Este documento registra, em ordem cronológica de descoberta, os 8 itens tratados:
1. reformulação do gerador de modelo para refletir a Figura 5.3;
2. bug no esquema gráfico (`system_schematic.png`);
3. decisão de design sobre `num_conductors`;
4. bug dimensional na matriz de impedância interna;
5. bug dimensional na matriz de retorno pelo solo (quasi-TEM);
6. `conductor_order` do MATLAB sem o ECC (dados descartados silenciosamente);
7. índices de ECC errados no `PLOT_CONFIG` para as matrizes de retorno pelo solo;
8. matrizes de retorno pelo solo do MATLAB recortadas em vez de reordenadas.

Mais uma pendência conhecida e não corrigida (item 9, `utils/comsol_data.py`).

---

## Item 1 — Reformulação de `flat_scc_with_ecc_cable_model` para a Figura 5.3

**Onde:** `models/single_core_cable.py:436`

**O quê:** o método era uma cópia de `hdpe_shared_enclosed_model` e assumia um duto HDPE
compartilhado entre o SCC e o ECC (com o ECC "entalado" entre a parede do duto e a
superfície do SCC via `_calculate_ecc_center_trig`). A Figura 5.3, porém, não tem duto —
os cabos são diretamente enterrados.

**Correção:** reescrito para posicionar os 3 cabos SCC em arranjo plano (mesma lógica de
`conventional_three_phase_flat`) e o ECC próximo ao terceiro cabo, **sem encostar nele**,
controlado por dois parâmetros de `arrangement`:
- `ecc_alignment`: `'center'` (mesmo eixo horizontal dos centros dos SCC, usado em
  `andreata_case3.json`) ou `'bottom_tangent'` (linha tangente ao fundo dos cabos SCC,
  simulando uma vala comum);
- `ecc_horizontal_gap`: folga entre as superfícies externas do 3º cabo SCC e do ECC
  (`andreata_case3.json` usa `0.050` m — valor de referência, sem cotação exata na
  figura original, ajustável).

Toda a lógica de duto/cálculo trigonométrico do "canto" foi removida deste método (ela
permanece intacta em `hdpe_shared_enclosed_model`, usada por outros casos).

---

## Item 2 — Marcador `'+'` fora de posição em `system_schematic.png`

**Onde:** `mtl_main/graphics.py`

**O quê:** o título do gráfico aparecia como o genérico "Multicondutor Transmission
Line" e o marcador `'+'`/linha de cota de profundidade ficavam colados à superfície do
solo, acima da posição real dos cabos.

**Causa raiz (duas camadas):**
1. `_calculate_schematic_parameters` (linha ~93-113) não reconhecia o tipo
   `'scc-flat-ecc'`, caindo no ramo genérico (`h_factor = 1`, título padrão).
2. Com `h_factor == 1`, `_plot_conductor_graphic` **pula** o reposicionamento
   esquemático (mantém os cabos na profundidade real), mas `_schematic_annotations`
   (linha 330) **sempre** aplicava o reposicionamento — os dois ficavam
   inconsistentes entre si sempre que um tipo caísse no ramo genérico.

**Correção:**
- `mtl_main/graphics.py:105` — adicionado `'scc-flat-ecc'` como tipo reconhecido
  (`h_factor = -2.5`, mesma convenção de `'scc'`; título "Flat-Buried SCC Cables with
  ECC").
- `mtl_main/graphics.py:330` (`_schematic_annotations`) — passou a espelhar o mesmo
  bypass de `h_factor == 1` que `_plot_conductor_graphic` já tinha, para que essa classe
  de inconsistência não volte a aparecer caso outro tipo futuro caia no ramo genérico.

**Validado:** `system_schematic.png` regenerado — cruz alinhada ao centro dos cabos, na
profundidade correta (h = 1,20 m).

---

## Item 3 — Decisão de design: `num_conductors` fixo em 3, fora do JSON

**Onde:** `models/single_core_cable.py:469`

**Discussão:** o método lia `num_scc_conductors = self.arrangement.get('num_conductors',
3)` do JSON. Ficou em aberto se esse número deveria contar só os cabos SCC (3) ou
incluir o ECC (4).

**Decisão:** `num_conductors` no JSON foi **removido**; `num_scc_conductors` passou a
ser um literal `3` no código, já que o método é inerentemente trifásico (mesma
convenção de `conventional_three_phase_flat`, que também não lê esse campo do JSON — a
cardinalidade está no próprio método, não em configuração externa). O ECC nunca é
contado por este campo em nenhum caso do repositório (ver `hdpe_ecc_2000mm2.json` /
`hdpe_ecc_300mm2.json`, onde `num_conductors: 1` refere-se só ao SCC hospedeiro).

---

## Item 4 — Matriz de impedância interna: `None + None` / dimensão incompatível

**Onde:** `analytical_forms/single_core_cable.py` (`InternalPerUnitParameters`) e
`mtl_main/strategy.py` (`SingleCoreCableWithECCStrategy`)

**Sintoma original:**
```
TypeError: unsupported operand type(s) for +: 'NoneType' and 'NoneType'
```
em `parameters_approximation`, linha `'Zcs': z11 + z12 + z2i`.

**Causa raiz (três camadas):**

1. **Formato errado de `context.scc`.** `SingleCoreCableWithECCStrategy._extract_scc_parameters`
   (`mtl_main/strategy.py:1056`) agrupa condutores por cabo físico e retorna um dict
   `{cp_key: {core_..., sheath_...}}` — correto para permitir cabos heterogêneos, mas
   incompatível com `InternalPerUnitParameters.parameters_by_bessel`/`parameters_approximation`,
   que esperavam `self.model.scc` como um dict **achatado** de uma única seção
   transversal (`scc['core_outer_radius']` etc. no nível raiz). Toda checagem
   `'core_outer_radius' in scc` dava `False` (as chaves eram strings de posição, não
   nomes de parâmetro) → `z11`, `z12`, `z2i` ficavam `None`.

2. **Montagem da matriz assumia N cabos idênticos.** `matrices()` (antiga versão)
   calculava um único bloco `Zij` (M×M) e montava a matriz completa via
   `np.kron(np.identity(N), Zij)` — não há como isso representar 3 cabos com M=2
   (núcleo+blindagem) e 1 cabo com M=1 (ECC) ao mesmo tempo.

3. **O ECC nem era capturado.** O laço de agrupamento de `_extract_scc_parameters`
   só aceitava `name in ('core', 'sheath', 'armor')` — o condutor `'ecc'` era
   descartado (`continue`), então seus dados geométricos nunca chegavam a
   `context.scc`.

**Correção:**
- `mtl_main/strategy.py:1056-1108` (`SingleCoreCableWithECCStrategy._extract_scc_parameters`)
  — passou a agrupar também condutores `'ecc'`; quando não há `core` numa posição (ou
  seja, é a posição do próprio ECC), o condutor ECC assume o papel de `core` nas
  fórmulas (fisicamente correto: um ECC isolado é só um condutor maciço com sua própria
  isolação, sem bainha/blindagem — mesma formulação do ramo "SCC com core,
  core_insulation").
- `analytical_forms/single_core_cable.py:832` (`_build_cable_block`) — lógica de
  montagem de `Zij_values`/`Pij` extraída para um método reutilizável por cabo
  (recebe `scc` explicitamente, em vez de sempre ler `self.model.scc`).
- `parameters_by_bessel`, `parameters_approximation`, `parameters_hybrid` — passaram a
  aceitar um parâmetro opcional `scc` (default: `self.model.scc`, preservando o
  comportamento anterior para os tipos homogêneos).
- `analytical_forms/single_core_cable.py:944` (`matrices`) — quando
  `self.model.mtl_type == 'scc-flat-ecc'`, delega para o novo método
  `_matrices_heterogeneous` (linha 989): itera cada grupo de cabo em `self.model.scc`
  (3 blocos 2×2 dos SCC + 1 bloco 1×1 do ECC) e monta a matriz final **diagonal em
  bloco** — sem `np.kron` uniforme —, já que não há acoplamento interno entre
  condutores de cabos físicos diferentes nesta etapa (o acoplamento entre cabos só
  entra depois, no retorno pelo solo).
- `analytical_forms/single_core_cable.py:467` (`_sum_or_none`) — corrigido de
  passagem um bug latente pré-existente, agora alcançável: o campo de conveniência
  `'Zcs'` fazia `z11 + z12 + z2i` incondicionalmente, mesmo quando `z2i` (que só existe
  se houver bainha) é `None` — exatamente o caso do ECC. Esse campo nunca é lido em
  nenhum lugar do pipeline; a soma agora ignora termos `None` em vez de lançar exceção.

**Validado:** matriz interna resultante tem shape `(num_freq, 7, 7)`, sem `NaN`, com os
blocos fora da diagonal entre cabos diferentes exatamente zero (desacoplamento interno
confirmado numericamente). Regressão limpa em `andreata_case1`, `scc_34kV_andreata`,
`hdpe_2000mm2`, `hdpe_300mm2`, `hdpe_ecc_2000mm2`, `scc_138kV_prysmian`, `scc_flat_xue`
(tipos homogêneos `scc`, `hdpe`, `shared-hdpe` não afetados).

---

## Item 5 — Matriz de retorno pelo solo (quasi-TEM): erro de broadcast 4×4 → 7×7

**Onde:** `analytical_forms/single_core_cable.py:1183`
(`PerUnitParameters.quasi_tem_approx_matrices`)

**Sintoma original:**
```
ValueError: could not broadcast input array from shape (4,4) into shape (7,7)
```
em `Zg[i, :, :] = np.kron(z0_jk[i, :, :], np.ones((M, M)))`.

**Causa raiz:** mesmo padrão do Item 4, uma camada acima. `z0_jk` (impedância de
retorno pelo solo) já vinha corretamente calculada como 4×4 — uma linha/coluna por
cabo físico (3 SCC + 1 ECC), pois `earth_return_parameters` usa `num_sc_cables`
(contagem por posição, sempre correta). O problema era só a expansão: `M =
self.model.num_conductors_per_scc` é uma média inteira (`total_conductores //
num_cabos` = 7 // 4 = 1) que assume todo cabo ter o mesmo número de condutores —
`np.kron(4×4, ones(1,1))` produz 4×4, incompatível com a `Zi` 7×7 do Item 4.

**Correção:**
- `analytical_forms/single_core_cable.py:478` — nova função `_expand_by_block_sizes(matrix,
  block_sizes)`: generaliza `np.kron(matrix, np.ones((M, M)))` para blocos de tamanho
  variável, via `np.repeat` (linhas e depois colunas) com uma lista de tamanhos por
  cabo em vez de um M uniforme. Quando todos os tamanhos são iguais a M, o resultado é
  idêntico ao `kron` original — ou seja, é uma generalização estrita, segura para os
  tipos homogêneos também.
- `InternalPerUnitParameters.matrices()` e `_matrices_heterogeneous()` passaram a
  incluir `'block_sizes'` no dicionário retornado (`[M]*N` no caminho homogêneo,
  `[2, 2, 2, 1]` no `scc-flat-ecc`) — uma única fonte de verdade sobre quantos
  condutores cada cabo contribui, na mesma ordem das matrizes de retorno pelo solo.
- `quasi_tem_approx_matrices` (linha 1193) passou a ler `internal_matrices['block_sizes']`
  em vez de `self.model.num_conductors_per_scc`, chamando `_expand_by_block_sizes` para
  montar `Zg`/`Pg` — sem mais laço por frequência, já que `np.repeat` opera direto no
  array 3D `(freq, N, N)`.

**Validado:** `Zg` resultante tem shape `(num_freq, 7, 7)`; conferido numericamente que
o bloco (cabo 1 × ECC) reproduz exatamente `z0_jk[:, 0, 3]` e o bloco (cabo 1 × cabo 2)
reproduz `z0_jk[:, 0, 1]`, ambos com diferença máxima `0.0`. `andreata_case3.py` roda
até o fim sem erro. Regressão limpa nos mesmos 7 casos homogêneos do Item 4.

---

## Item 6 — `conductor_order` do MATLAB sem o ECC (dados descartados silenciosamente)

**Onde:** `testData/andreata_case3/andreata_case3.py` (chamada a `MatlabDataReader.get_scc_scenario_data`)

**Sintoma original:** três warnings ao rodar, sem crash:
```
Warning: MATLAB data not found for 'measured' with key 'series_impedance_matrix'.
Warning: MATLAB data not found for 'measured' with key 'internal_impedance_matrix'.
Warning: MATLAB data not found for 'measured' with key 'internal_admittance_matrix'.
```

**Causa raiz:** os arquivos `.mat` de referência já vinham no formato cheio 7×7 (confirmado
via `scipy.io.loadmat`: shape `(7, 7, 90)`), mas `conductor_order=[0, 3, 1, 4, 2, 5]`
— herdado do `andreata_case1`, que não tem ECC — só listava 6 índices.
`MatlabDataReader._reorder_conductor_matrix` (`utils/matlab_data.py`) faz
`matrix[:, conductor_order, :][:, :, conductor_order]`: com apenas 6 índices, a
linha/coluna do ECC (índice 6) era silenciosamente descartada, produzindo uma matriz
`'measured'` 6×6. Os componentes de plotagem que pedem `p=6, q=6` (termo próprio do ECC)
então não encontravam esse índice (`IndexError`, capturado pelo `except` genérico do
plotter — `scc_plotter.py:69` — e impresso apenas como "not found").

**Correção:** `conductor_order=[0, 3, 1, 4, 2, 5, 6]` — o ECC não tem par núcleo/bainha
para trocar de posição, então permanece no fim (índice 6) tanto na convenção MATLAB
(tipo-agrupada) quanto na convenção pyLCP (cabo-agrupada).

**Validado:** os três warnings desaparecem; `andreata_case3.py` roda sem alterar as
demais curvas (índices 0-5 inalterados pela mudança).

---

## Item 7 — Índices de ECC errados no `PLOT_CONFIG` para as matrizes de retorno pelo solo

**Onde:** `testData/andreata_case3/plot_config.py`

**Sintoma original:**
```
IndexError: index 6 is out of bounds for axis 1 with size 4
```
em `scc_plotter.py:32` (`_plot_scc_matrix`), ao plotar `earth_return_impedance_ecc`.

**Causa raiz:** as configs `earth_return_impedance_ecc`, `earth_return_admittance_ecc` e
`earth_return_potential_coeff_ecc` usavam `path: ['earth_return_parameters', ...]` com
`p=6, q=6` — mas `PerUnitParameters.earth_return_parameters()`
(`analytical_forms/single_core_cable.py:1081`) retorna matrizes indexadas **por cabo
físico** (`N = num_sc_cables = 4`: fase A, fase B, fase C, ECC — ver Item 5), não por
condutor. Nesse espaço o ECC é o índice 3, não 6. Já as configs
`mutual_impedance_phase_a_sheath_ecc`, `self_impedance_ecc` e `self_admittance_ecc`, que
usam `path: ['quasi_tem_matrices', ...]`, já estavam corretas com `p=6, q=6`, pois essas
matrizes foram expandidas para o espaço por-condutor (N=7) via `_expand_by_block_sizes`
(Item 5).

**Correção:**
- `earth_return_admittance_ecc`: mantido em `path: ['earth_return_parameters',
  'admittance_matrix']`, com `p=3, q=3` (índice do ECC no espaço por-cabo). Não há
  equivalente por-condutor válido para essa matriz — ver observação abaixo.
- `earth_return_impedance_ecc` / `earth_return_potential_coeff_ecc`: migrados para
  `path: ['quasi_tem_matrices', 'earth_return_impedance_matrix']` /
  `['quasi_tem_matrices', 'earth_return_potential_coefficient']`, com `p=6, q=6` — essas
  chaves já existem em `quasi_tem_approx_matrices` (linha 1219-1226), expandidas para
  N=7, e batem índice-a-índice com o MATLAB reordenado (ver Item 8).

**Observação (não corrigida):** `quasi_tem_approx_matrices` retorna
`'earth_return_admittance_matrix': Yg`, mas `Yg` nunca é preenchido (fica
`np.zeros_like(Zi, dtype=complex)`, dead code pré-existente) — por isso
`earth_return_admittance_ecc` não pôde ser migrado como as outras duas. Como essa config
não tem `matlab_series_to_plot`, isso não gera warning nem crash hoje; fica registrado
para quando alguém for calcular `Yg` de fato.

**Validado:** `andreata_case3.py` roda até o fim sem `IndexError`.

---

## Item 8 — Matrizes de retorno pelo solo do MATLAB recortadas em vez de reordenadas (ECC descartado)

**Onde:** `utils/matlab_data.py` (`MatlabDataReader.get_scc_scenario_data`)

**Sintoma original (após o Item 7):**
```
Warning: MATLAB data not found for 'measured' with key 'earth_return_impedance_matrix'.
Warning: MATLAB data not found for 'measured' with key 'earth_return_potential_coefficient_matrix'.
```

**Causa raiz:** `earth_return_impedance_matrix` e `earth_return_potential_coefficient_matrix`
recebiam um tratamento diferente das outras 4 matrizes: em vez de
`_reorder_conductor_matrix(..., conductor_order)`, o código recortava apenas o bloco
superior-esquerdo `[:num_phases, :num_phases]` (3×3), sob a suposição — correta para
`andreata_case1` (sem ECC), mas desatualizada para `andreata_case3` — de que essas
matrizes só existem "a nível de cabo/fase". Confirmado via `scipy.io.loadmat` que os
`.mat` de retorno pelo solo já vêm no mesmo formato cheio 7×7, tipo-agrupado-com-redundância,
das demais 4 matrizes: `M[0,0]==M[3,3]` (core_A==sheath_A), `M[1,1]==M[4,4]`,
`M[2,2]==M[5,5]`, e `M[6,6]` é o valor próprio do ECC. O recorte a 3×3 descartava esse
índice 6 por completo, então `'measured'` nunca tinha dado para o ECC.

**Correção:** as duas matrizes passaram a usar `_reorder_conductor_matrix(matrix,
conductor_order)`, igual às outras 4 — sem recorte especial. O parâmetro `num_phases`
(que só servia a esse recorte) foi removido de `get_scc_scenario_data`, e as chamadas em
`andreata_case1.py`/`andreata_case3.py` atualizadas.

**Validado:**
- Os dois warnings desaparecem; nenhum warning novo surge em `andreata_case1` (os índices
  `p=0, q=0` de fase-A apontam para o mesmo valor de `core_A` tanto no formato recortado
  quanto no formato cheio reordenado — confirmado numericamente).
- Checagem numérica pontual (índice de frequência 45, ≈2,15 kHz): `Zg` pyLCP
  `[:,6,6] = 3,532e-4 + j4,951e-3` vs. MATLAB `3,531e-4 + j4,952e-3`; `Pg` pyLCP
  `= 5,580e4 + j5,357e5` vs. MATLAB `5,580e4 + j5,354e5` — mesma ordem de grandeza,
  consistente com dado medido vs. formulação analítica (mesmo padrão de `Zi_77`/`Yi_77`,
  Item 4).

---

## Item 9 (pendente) — Mesmo padrão de bug em `utils/comsol_data.py`

**Onde:** `utils/comsol_data.py:731` (`ComsolPostProcessor.get_quasi_tem_approx_matrices`)

**O quê:** função irmã de `quasi_tem_approx_matrices` (usada para comparar com dados
COMSOL), com o mesmo padrão de expansão por `np.kron(z0_jk, np.ones((M, M)))` — porém
com `M = 2` **hardcoded** (linha 738), nem sequer derivado do modelo.

**Por que não quebrou ainda:** `andreata_case3.py` só chama esse método dentro do bloco
`if cmsl_params is not None:` — e não há arquivo `Results/cmsl_ground_return_impedance.txt`
para este caso ("Processamento COMSOL ignorado (dados não disponíveis)" no console), então
o caminho nunca foi exercitado.

**Risco:** se um dado COMSOL for adicionado a este caso no futuro, vai quebrar do mesmo
jeito que os Itens 4 e 5 quebraram — `np.kron(4×4, ones((2,2)))` dá 8×8, incompatível
com a `Zi` 7×7.

**Correção proposta (não aplicada):** como `ComsolPostProcessor` não tem acesso a
`self.model`, a correção exigiria passar `block_sizes` como parâmetro explícito à
função (em vez de derivá-lo de `self.model.num_conductors_per_scc` como hoje), reusando
a mesma `_expand_by_block_sizes` do Item 5. Fica registrado aqui para quando houver dado
COMSOL disponível para este caso.
