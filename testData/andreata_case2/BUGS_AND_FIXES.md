# Diagnósticos e correções — `andreata_case2` (Configuração 2, Figura 5.2)

**Contexto:** `andreata_case2` modela a Configuração 2 da referência (Figura 5.2): três
cabos SCC de potência (núcleo + blindagem) em arranjo plano, cada um instalado dentro do
seu próprio duto de HDPE, enterrados no solo. Diferente de `andreata_case1`
(Configuração 1, cabos diretamente enterrados) e `andreata_case3` (Configuração 3, SCC +
ECC sem duto), aqui **não existe solução analítica de retorno à terra que modele o
duto** — a formulação de Zg/Yg do pyLCP enxerga só a posição/raio externo de cada cabo,
nunca a presença de um duto. A validação do efeito do duto depende inteiramente de
COMSOL/MATLAB.

O script encontrado no início desta rodada de trabalho era de uma **origem de
desenvolvimento mais antiga**: modelava um único cabo dentro de um duto (sem as três
fases, sem acoplamento entre cabos, sem retorno à terra), usando uma API do COMSOL
própria de cabo isolado e um plotter (`HDPEPlotter`) e `plot_config.py` num formato
diferente do adotado por `andreata_case1`/`andreata_case3`. Essa mesma linhagem antiga
ainda serve os casos `hdpe_2000mm2`, `hdpe_300mm2`, `hdpe_ecc_2000mm2`, `hdpe_ecc_300mm2`
e não foi alterada.

Este documento registra, em ordem cronológica, o que foi feito:
1. generalização do gerador de modelo para 3 fases com duto individual (Figura 5.2);
2. bug de posicionamento — duto concêntrico em vez de excêntrico;
3. bug de aliasing entre os dutos das 3 fases;
4. atualização do JSON para o arranjo trifásico;
5. reescrita de `andreata_case2.py` no molde de `andreata_case1`/`andreata_case3`;
6. generalização dos índices fixos de blindagem (`model_2[2]`/`model_3[2]`) para N fases;
7. reescrita de `plot_config.py` no formato dirigido a config do `SCCPlotter`;
8. validação cruzada contra o MATLAB do `andreata_case1` (cenário sem duto);
9. limitações conhecidas e pendências de dados externos.

---

## Item 1 — Novo gerador de modelo: `flat_hdpe_enclosed_model`

**Onde:** `models/single_core_cable.py`

**O quê:** não existia nenhum método que gerasse N cabos em arranjo plano, cada um
dentro do seu próprio duto de HDPE — `eccentric_hdpe_enclosed_model`/
`concentric_hdpe_enclosed_model` eram hardcoded para 1 cabo só (a origem do script
antigo). A camada de MTL (`SingleCoreCableInHDPEStrategy` em `mtl_main/strategy.py`), em
contraste, já era genérica para N cabos — ela é praticamente uma cópia de
`SingleCoreCableStrategy` (usada pelo `type: "scc"` trifásico do case1), agrupando
condutores por `center_point`. Ou seja, a limitação estava só na geração de geometria,
não na infraestrutura de cálculo.

**Correção:** novo método `flat_hdpe_enclosed_model(host_conductor='sheath')`,
generalizando `eccentric_hdpe_enclosed_model` (1 cabo) para N cabos, do mesmo modo que
`underground_flat_model` generaliza `conventional_single_phase`: mesmos
`burial_depth`/`spacing`/`num_conductors` do `arrangement`, um duto idêntico injetado em
cada cabo.

---

## Item 2 — Duto gerado concêntrico; deveria ser excêntrico (cabo no fundo do duto)

**Onde:** `models/single_core_cable.py` (`flat_hdpe_enclosed_model`)

**O quê:** a primeira versão do método usava geometria **concêntrica** (cabo e duto
compartilhando o mesmo centro) — por analogia direta com `concentric_hdpe_enclosed_model`
e com uma leitura inicial da Figura 5.2. O usuário apontou que essa não é a convenção
usada no restante do repositório: `schematic_original_hdpe.png` (caso `hdpe_300mm2`, já
existente) mostra o cabo **apoiado na superfície interna inferior do duto** — a mesma
convenção de `eccentric_hdpe_enclosed_model` —, e não centrado nele.

**Impacto:** com `burial_depth` (`h = 1,20 m`) fixado no centro do **cabo** (não do
duto), a geometria concêntrica colocava o centro do duto exatamente na mesma cota do
cabo; a correta exige deslocar o centro do duto **para cima**, preservando `h` no centro
do cabo.

**Correção:** reescrito para excentricidade, replicando a lógica de
`eccentric_hdpe_enclosed_model`: `vertical_offset = enclosure_inner_radius -
cable_outer_radius`; o centro do duto de cada fase é `(cable_center[0], cable_center[1] +
vertical_offset)`. Aplicado identicamente às 3 fases (mesmo duto, mesmo deslocamento).

**Validado:** `system_schematic.png` regenerado — os 3 cabos aparecem apoiados no fundo
de seus respectivos dutos, batendo visualmente com `schematic_original_hdpe.png` e com a
Figura 5.2 (cabos em x = 0 / 0,2 / 0,4 m, profundidade do centro do cabo h = 1,20 m).

---

## Item 3 — Bug de aliasing: os 3 dutos compartilhariam o mesmo `center_point`

**Onde:** `models/single_core_cable.py` (`flat_hdpe_enclosed_model`)

**O quê:** o padrão usado pelos métodos de 1 cabo (`eccentric_/concentric_hdpe_enclosed_model`)
injeta o `enclosure_data` **por referência** direta no condutor host encontrado (só há
um). Generalizado ingenuamente para N cabos, o mesmo dict `enclosure_data` (obtido uma
única vez de `self.sheath['enclosure']`) seria compartilhado pelas 3 blindagens — cada
atribuição de `v['enclosure']['center_point']` no laço sobrescreveria a mesma instância,
de forma que **todas as 3 fases acabariam com o `center_point` do duto da última fase
processada**.

**Correção:** cada condutor recebe sua própria cópia via `copy.deepcopy(enclosure_data)`
antes de setar seu `center_point` individual. (`import copy` adicionado ao topo do
arquivo.)

---

## Item 4 — `andreata_case2.json`: arranjo trifásico

**Onde:** `testData/andreata_case2/andreata_case2.json`

**O quê:** `arrangement` descrevia 1 cabo (`"type": "single"`, `"num_conductors": 1`,
`"spacing": null`, `"fourier_order": 4`).

**Correção:** atualizado para `"type": "flat"`, `"num_conductors": 3`, `"spacing": 0.2`
(m, centro-a-centro, conforme Figura 5.2), `"fourier_order": 10` (alinhado ao
`andreata_case1`, antes divergente sem motivo aparente). `burial_depth` (1,2 m) e a
definição de `cable_definition` (cabo + duto) já correspondiam à Figura 5.2 e não foram
alteradas. `name`/`note` atualizados para descrever a Configuração 2 trifásica.

---

## Item 5 — Reescrita de `andreata_case2.py` no molde do case1/case3

**Onde:** `testData/andreata_case2/andreata_case2.py`

**O quê:** o script original refletia a linhagem antiga de desenvolvimento (ver seção de
Contexto) e não podia simplesmente ser "adaptado" — várias peças inteiras estavam
ausentes porque, com 1 cabo só, não faziam sentido:
- nunca chamava `apply_semiconducting_layer_correction`, apesar do JSON já definir
  `semiconducting_layer` para núcleo e blindagem (correção nunca aplicada — inconsistência
  latente, independente da mudança para 3 fases);
- não calculava retorno à terra (`PerUnitParameters`) nem matrizes quasi-TEM — sem
  acoplamento entre fases, isso nunca tinha sido necessário;
- não usava `MatlabDataReader` (case1/case3 usam);
- lia dados COMSOL pela API de cabo isolado (`cmsl_coaxial_cable_impedance`/
  `get_coaxial_cable_parameters`/`get_scc_internal_impedance_elements`), incompatível
  com a família de retorno à terra trifásica (`cmsl_ground_return_impedance`/
  `get_earth_return_parameters`) usada por case1/case3;
- plotava com `HDPEPlotter` (API de "data_series" própria, incompatível com o formato de
  `plot_config.py` dirigido a config usado por `SCCPlotter`).

**Correção:** reescrito seguindo o esqueleto do `andreata_case1.py`/`andreata_case3.py`:
- `model_0 = flat_hdpe_enclosed_model()` → `apply_semiconducting_layer_correction(...)` →
  `MulticonductorTransmissionLine(model_0)` — geometria física real, usada para o
  esquemático e como base do ERS/GMD (via `EquivalentRadiiSystems(mtl_0)`, que fornece
  `r4`..`r7`);
- `model_1 = underground_flat_model()` (+ mesma correção semicondutora) — cenário
  "Underground", ignora o duto por completo; serve de baseline e de contraparte direta
  do `andreata_case1` (ver Item 8);
- `model_2`/`model_3` — os dois modelos equivalentes de duto já existentes no script
  antigo (ERS ponderado por área; GMD caso 3.1) foram preservados como o diferencial
  físico do caso, agora extraídos para uma função compartilhada `_override_sheath_insulation`
  (ver Item 6) e alimentando os mesmos 3 cenários (`'1'`, `'2'`, `'3'`) usados no resto do
  pipeline;
- bloco de retorno à terra adicionado ao laço de cenários: `PerUnitParameters(...).earth_return_parameters(...)`
  + `quasi_tem_approx_matrices(...)`, análogo ao case1, mas com uma única formulação de
  solo (`magalhaes_xue`) — o eixo de comparação relevante aqui é bare/ERS/GMD, não
  formulações de solo (esse é o assunto do case1);
- `MatlabDataReader` adicionado (prefixo `andreata_hdpe`, ver Item 9 sobre a suposição de
  nome de arquivo);
- COMSOL trocado para a família de retorno à terra (`cmsl_ground_return_impedance` /
  `get_earth_return_parameters('rho_g_100_epsr1_1_mf')` — chave escolhida por já
  corresponder ao solo do JSON, ρ=100 Ω·m / εr=1); os parâmetros internos usados para
  compor a matriz quasi-TEM do COMSOL vêm do modelo GMD caso 3.1 (`mtl_3`), por ser a
  variante mais completa das três — só o retorno à terra é, de fato, medido no COMSOL;
- `HDPEPlotter` trocado por `SCCPlotter` (`compare_complete_matrices`), com
  `plot_config.py` reescrito (ver Item 7);
- `compare_capacitance()` (tabela impressa de C₁₁/C₂₂ da fase A) preservada quase
  inalterada — ela já indexava por cenário (`'1'`/`'2'`/`'3'`), então generaliza sem
  mudanças para N=3 fases (índices `[0,0]`/`[1,1]` da matriz de capacitância continuam
  sendo núcleo/blindagem da fase A, já que a ordenação por cabo não muda);
- leitura de `get_shunt_capacitance_elements()` (COMSOL) agora gateada explicitamente por
  `'cmsl_shunt_params' in cmsl_processor.cmsl_reader.data` — no script antigo, essa
  chamada (que acessa `self.cmsl_reader.data['cmsl_shunt_params']` diretamente, sem
  tratamento de ausência) era gateada pela presença de um arquivo COMSOL *diferente*
  (`cmsl_coaxial_cable_impedance`), o que já era impreciso e teria virado um `KeyError`
  não tratado assim que a família COMSOL fosse trocada.

---

## Item 6 — Índices fixos de blindagem (`model_2[2]`) não generalizavam para N fases

**Onde:** `testData/andreata_case2/andreata_case2.py`

**O quê:** o script antigo (1 cabo) localizava a blindagem para aplicar a permissividade
equivalente por posição fixa: `model_2[2]['enclosure'] = None`,
`model_2[2]['insulation']['thickness'] = ...` (índice `2` = única blindagem existente,
já que a ordem de inserção é `0=solo, 1=núcleo, 2=blindagem`). Com 3 fases, as blindagens
ficam nos índices `2, 4, 6` — o código antigo, se simplesmente reaproveitado, teria
alterado só a blindagem da fase A e deixado B/C com o duto físico intacto (inconsistência
silenciosa, sem erro).

**Correção:** extraído para a função `_override_sheath_insulation(model, thickness,
relative_permittivity)`, que varre `model.values()` por `conductor_name == 'sheath'` e
aplica a substituição a todas as blindagens encontradas — correto para qualquer N.
`model_3` mantém a espessura original da isolação (`model_0[2]['insulation']['thickness']`,
lida antes da substituição — Caso 3.1 mantém r0 = r5), alterando só a permissividade;
`model_2` usa a espessura estendida até r7.

---

## Item 7 — Reescrita de `plot_config.py` no formato do `SCCPlotter`

**Onde:** `testData/andreata_case2/plot_config.py`

**O quê:** o arquivo antigo usava um esquema `data_series` (lista de dicts com
`source`/`scenario_key`/`data_path`/`plot_style`/`series`) próprio do `HDPEPlotter`,
incompatível com o formato `series_to_plot`/`path`/`p`/`q`/`components` que o `SCCPlotter`
espera (mesmo usado por `andreata_case1`/`andreata_case3`).

**Correção:** reescrito nesse segundo formato. Templates de série:
`DUCT_MODEL_TEMPLATE` (cenários `'1'`/`'2'`/`'3'` — Underground/ERS/GMD, uma linha cada),
`COMSOL_TEMPLATE` (`rho_g_100_epsr1_1_mf`), `MATLAB_TEMPLATE` (`'measured'`, dado
específico do caso, ainda não existente). 9 gráficos via `compare_complete_matrices`:
6 internos (núcleo/blindagem/mútua × impedância/admitância, todos sensíveis ao modelo de
duto porque ERS/GMD alteram espessura e/ou permissividade da isolação externa da
blindagem) + 3 de retorno à terra por fase (impedância, admitância, coeficiente de
potencial — ver nota no Item 5 sobre a formulação não modelar o duto).

---

## Item 8 — Validação cruzada com o MATLAB do `andreata_case1`

**Onde:** `testData/andreata_case2/andreata_case2.py` (`validate_against_case1_reference`),
`testData/andreata_case2/plot_config.py` (`MATLAB_CASE1_NO_DUCT_TEMPLATE`)

**Motivação:** o cenário `'1'` ("Underground", duto ignorado, gerado por
`underground_flat_model()`) descreve exatamente a mesma geometria física do
`andreata_case1` (3 cabos SCC em arranjo plano, sem duto). Como o `andreata_case1` já tem
dados MATLAB de referência validados (ver `andreata_case1/BUGS_MATLAB_IMPORT.md`), eles
servem de checagem cruzada de código para a parte do pipeline que os dois casos
compartilham (`InternalPerUnitParameters`, `PerUnitParameters`), sem depender dos dados
COMSOL/MATLAB específicos do duto (que ainda não existem).

**Implementação:** um segundo `MatlabDataReader`, apontado para o caminho de script do
`andreata_case1` (não o `__file__` do próprio case2) — assim ele lê
`testData/andreata_case1/Results/*.mat` em vez de `testData/andreata_case2/Results/*.mat`
— carrega o mesmo prefixo (`'andreata'`) e `conductor_order` (`[0,3,1,4,2,5]`) que o
`andreata_case1.py` usa para si mesmo. O resultado é armazenado como o cenário
`'case1_no_duct'` dentro de `pul_data['matlab']['scenarios']`, ao lado do cenário
`'measured'` (dados próprios do case2, com duto, ainda ausentes) — os dois podem
coexistir e ser plotados juntos sem conflito, já que o `SCCPlotter` indexa por chave de
cenário.

`validate_against_case1_reference()` imprime o erro relativo máximo (norma-∞) entre a
matriz quasi-TEM completa do cenário `'1'` (`series_impedance_matrix`,
`shunt_admittance_matrix`) e a referência MATLAB do case1.

**Resultado (com os dados MATLAB do `andreata_case1` já existentes em seu `Results/`):**
```
Zs  (impedância série completa):  erro relativo máx. (norma inf) = 0.06%
Ysh (admitância shunt completa):  erro relativo máx. (norma inf) = 0.15%
```
Concordância excelente — confirma que a geração do modelo bare trifásico e o cálculo de
parâmetros internos/retorno à terra estão corretos nesta parte do pipeline. **Não** valida
o efeito do duto em si (`'2'`/`'3'`), que não tem contraparte no `andreata_case1`.

`plot_config.py` ganhou `MATLAB_CASE1_NO_DUCT_TEMPLATE` (marcador `+` laranja, chave
`'case1_no_duct'`), somado a `MATLAB_TEMPLATE` nos gráficos internos e de retorno à terra,
para sobrepor essa referência visualmente à linha preta do cenário `'1'`.

---

## Item 9 — Limitações conhecidas e pendências de dados externos

- **COMSOL do duto ainda não existe.** `cmsl_processor.get_general_parameters('cmsl_ground_return_impedance')`
  retorna `None` (arquivo ausente em `Results/`) — o bloco inteiro é ignorado
  graciosamente (`pul_data['comsol'] = {}`), igual ao `andreata_case1`. Quando adicionado,
  deve seguir o mesmo formato de colunas do `andreata_case1`
  (`cmsl_ground_return_impedance.txt`, sufixos `_rho_g_100_epsr1_1_mf_vcoil_1/2/3`) — a
  chave de cenário `'rho_g_100_epsr1_1_mf'` já está fixada no código combinando com o
  solo do JSON (ρ=100 Ω·m, εr=1); se o arquivo real usar outro `rho_g`/`epsr1`, ajustar a
  chave em `andreata_case2.py` e `COMSOL_TEMPLATE` (`plot_config.py`).
- **MATLAB do duto ainda não existe.** `MatlabDataReader(__file__, ...)` procura por
  `.mat` em `testData/andreata_case2/Results/` com prefixo `'andreata_hdpe'` (ex.:
  `andreata_hdpe_frequency_range.mat`, `andreata_hdpe_series_impedance_matrix.mat`, etc.
  — mesmo conjunto de variáveis que o `andreata_case1` usa, só com prefixo diferente).
  Esse prefixo foi escolhido para não colidir com o `'andreata'` puro reservado ao
  case1 (ver Item 8) — se o prefixo real vier diferente, ajustar a string em
  `andreata_case2.py`.
- **`cmsl_shunt_params.txt` existente em `Results/` é do estudo antigo de cabo único**,
  reaproveitado apenas para a linha "COMSOL (energy method)" da tabela de capacitância —
  não é uma referência do sistema trifásico e deve ser substituído quando os dados novos
  chegarem (mesmo nome de arquivo esperado, `get_shunt_capacitance_elements()` não muda).
  Os demais arquivos legados em `Results/` (`cmsl_coaxial_cable_impedance.txt`,
  `cmsl_series_impedance_core_excitation.txt`, `cmsl_series_impedance_sheath_excitation.txt`)
  não são mais lidos pelo script atual (API de cabo isolado, abandonada — ver Item 5);
  ficam como dado morto em `Results/` até uma limpeza futura, sem bloquear a execução.
- **Retorno à terra nunca modela o duto**, em nenhum dos 3 cenários — ver nota extensa no
  código (`andreata_case2.py`, dentro do laço de cenários) e no Item de Contexto acima.
  A única forma de validar o efeito real do duto no retorno à terra é a comparação
  COMSOL/MATLAB pendente.

---

## Módulos legados não alterados

`models/hdpe.py`, `models/model_generator.py`, `models/isolated_conductors.py` e a classe
`HDPEPlotter` (`plotter/scc_plotter.py`) continuam servindo os casos `hdpe_2000mm2`,
`hdpe_300mm2`, `hdpe_ecc_2000mm2`, `hdpe_ecc_300mm2` — são a linhagem antiga de fato, e
`andreata_case2` deixou de depender deles (usa exclusivamente
`models.single_core_cable.SingleCoreCableModelGenerator`, `analytical_forms.single_core_cable`
e `plotter.scc_plotter.SCCPlotter`, a mesma base do `andreata_case1`/`andreata_case3`).
