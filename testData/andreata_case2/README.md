# andreata_case2

## Propósito

Script de simulação e validação para a **Configuração 2** (Figura 5.2) da referência:
três cabos coaxiais de núcleo único (SCC — núcleo + blindagem) em arranjo plano,
**cada um instalado dentro do seu próprio duto de HDPE**, enterrados no solo. Calcula
os parâmetros por unidade de comprimento (PUL) — internos (núcleo/blindagem) e de
retorno à terra — comparando três formas analíticas distintas de representar o efeito
dielétrico do duto, e faz validação cruzada com os dados MATLAB do `andreata_case1`.

Diferente do `andreata_case1` (Configuração 1, cabos diretamente enterrados) e do
`andreata_case3` (Configuração 3, SCC + ECC sem duto), este caso lida com uma limitação
estrutural do pyLCP: **não existe solução analítica de retorno à terra que modele o
duto** — a formulação de `Zg`/`Yg` enxerga apenas a posição e o raio externo de cada
cabo, nunca a presença de um duto de HDPE ao redor dele. A validação de fato do efeito
do duto no retorno à terra depende inteiramente de dados COMSOL/MATLAB específicos do
duto, que **ainda não existem** para este caso (ver [Limitações e Dados
Pendentes](#limitações-e-dados-pendentes)).

Referência: Andreata, Luis Eduardo Batista. *Análise das Características de Propagação
e de Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não
Metálicos no Contexto de Parques Eólicos.* Programa de Pós-Graduação em Engenharia
Elétrica, UFMG, 2025. https://hdl.handle.net/1843/2086

---

## Configuração do Sistema (JSON)

Parâmetros carregados de `andreata_case2.json`:

| Parâmetro | Valor |
|---|---|
| Arranjo | Plano (`flat`), 3 cabos |
| Profundidade de enterramento | 1,2 m (centro do cabo) |
| Espaçamento entre cabos | 0,2 m (centro a centro) |
| Solo — condutividade | 0,01 S/m (ρ = 100 Ω·m) |
| Solo — permissividade relativa | 1,0 |
| Núcleo — raio externo | 10,325 mm — condutividade: 38 MS/m |
| Isolação do núcleo | XLPE, 10,675 mm, εr = 2,2 |
| Bainha — raios int./ext. | 21,0 / 21,8 mm — condutividade: 60 MS/m |
| Isolação da bainha | PVC, 2,2 mm, εr = 2,8 |
| Duto (HDPE) — raios int./ext. | 50,8 / 57,8 mm — εr = 2,35, excêntrico (cabo apoiado no fundo) |
| Ordem de Fourier | 10 |

O duto é definido em `cable_definition.sheath.enclosure` — cada uma das 3 fases recebe
um duto idêntico, gerado de forma excêntrica (`vertical_offset = enclosure_inner_radius
- cable_outer_radius`) para que o cabo fique apoiado na parede interna inferior do duto,
com o centro do cabo (não do duto) na profundidade `burial_depth`.

---

## Estrutura do Script

### Função auxiliar `_override_sheath_insulation` (linhas 28–42)

Substitui o duto físico (`enclosure`) de cada blindagem por uma camada de isolação
homogênea equivalente sobre o próprio condutor `sheath` — a técnica ERS/GMD usada para
incorporar o efeito dielétrico do duto de HDPE nas fórmulas padrão de parâmetros
internos núcleo+blindagem, que não têm noção própria de uma camada de duto. Varre
`model.values()` procurando `conductor_name == 'sheath'` e aplica a substituição a
todas as blindagens encontradas — generaliza para qualquer número de fases (nos índices
fixos `2, 4, 6` do modelo de 3 cabos).

### Função `compare_capacitance` (linhas 44–70)

Imprime uma tabela comparando `C₁₁` (núcleo) e `C₂₂` (blindagem) da fase A, em µF/km,
para os três modelos de duto (`'1'` Underground, `'2'` ERS, `'3'` GMD caso 3.1) e,
quando disponível, a referência COMSOL (método de energia, `cmsl_shunt_params.txt`).
Os índices `[0,0]`/`[1,1]` da matriz de capacitância correspondem sempre a
núcleo/blindagem da fase A, independentemente do número de fases (ordenação por cabo:
`np.kron(identity(N), Zij)`).

### Função `validate_against_case1_reference` (linhas 72–102)

Confere se o cenário `'1'` ("Underground", duto ignorado) concorda com a referência
MATLAB do `andreata_case1` — mesma geometria física de 3 cabos SCC em arranjo plano,
sem duto. Não valida o efeito do duto em si (sem contraparte no case1); serve como
checagem cruzada de código para a parte do pipeline compartilhada entre os dois casos
(`underground_flat_model`, `InternalPerUnitParameters`, `PerUnitParameters`). Calcula o
erro relativo máximo (norma-∞) entre `Zs`/`Ysh` quasi-TEM completos do cenário `'1'` e a
matriz MATLAB do case1.

### Função `main()` (linhas 104–248)

#### 1. Construção dos modelos (linhas 107–140)

Quatro modelos MTL são construídos a partir de `SingleCoreCableModelGenerator`:

| Modelo | Geometria | Papel |
|---|---|---|
| `model_0` / `mtl_0` | `flat_hdpe_enclosed_model()` — geometria física real, com duto | Base do esquemático e do ERS/GMD (fornece `r4`..`r7` via `EquivalentRadiiSystems`). Não é usado diretamente como cenário. |
| `model_1` / `mtl_1` | `underground_flat_model()` — ignora o duto por completo | Cenário `'1'` — baseline "Underground", contraparte direta do `andreata_case1` |
| `model_2` / `mtl_2` | isolação equivalente ponderada por área (ERS) via `_override_sheath_insulation` | Cenário `'2'` — "Area-weighted ERS" |
| `model_3` / `mtl_3` | isolação equivalente por GMD, caso 3.1 (Lafaia, 2015), mantendo `r0 = r5` | Cenário `'3'` — "GMD case 3.1" |

`EquivalentRadiiSystems(mtl_0)` calcula os raios equivalentes (`equiv_rel_permittivity_epsr_area_weighted`,
`equivalent_parameters_from_gmd`) usados para parametrizar `model_2` e `model_3`.

#### 2. Dicionário `pul_data` (linhas 142–154)

```
pul_data
├── 'frequencies'  → 90 pontos log-espaçados em [10⁻², 10⁷] Hz
├── 'comsol'       → {'scenarios': {'rho_g_100_epsr1_1_mf': {}}}  (preenchido se arquivo disponível)
└── 'scenarios'    → '1' (Underground), '2' (ERS), '3' (GMD caso 3.1)
    (mesma formulação de solo para os três: zg_form = yg_form = 'magalhaes_xue' —
    o eixo de comparação aqui é o modelo de duto, não a formulação de solo, que é
    o assunto do andreata_case1)
```

#### 3. Dados MATLAB (linhas 156–181)

Dois `MatlabDataReader` distintos:

- **Dados próprios do duto** (`prefix='andreata_case2'`) — aponta para
  `Results/*.mat` do próprio `andreata_case2` (arquivos nomeados
  `andreata_case2_*.mat`); armazenados como cenário `'measured'`.
- **Dados do `andreata_case1`** (`prefix='andreata_case1'`) — aponta explicitamente
  para o script `testData/andreata_case1/andreata_case1.py`, lendo
  `andreata_case1/Results/andreata_case1_*.mat` (já existentes e validados).
  Armazenados como cenário `'case1_no_duct'`, ao lado de `'measured'`, para checagem
  cruzada com o cenário `'1'`.

Cada caso (`andreata_case1`, `andreata_case2`, `andreata_case3`) usa seu próprio
prefixo de arquivo (`andreata_case{N}_*.mat`) para evitar ambiguidade entre os
diretórios `Results/` dos três casos — antes todos usavam o prefixo genérico
`andreata_*`, o que causava colisões de nome ao copiar/gerar dados de um caso para
outro.

Os dois `MatlabDataReader` usam a mesma permutação de condutores
`conductor_order=[0, 3, 1, 4, 2, 5]` (agrupamento por tipo no MATLAB → agrupamento por
cabo no pyLCP), herdada do `andreata_case1`.

#### 4. Processamento COMSOL — opcional (linhas 183–203)

Tenta carregar `Results/cmsl_ground_return_impedance.txt` (família de retorno à terra,
`get_earth_return_parameters`) via `ComsolPostProcessor`. **Ainda não existe** para este
caso (ver [Limitações](#limitações-e-dados-pendentes)) — o bloco é ignorado com aviso e
`pul_data['comsol']` esvaziado. Quando presente, os parâmetros internos usados para
compor a matriz quasi-TEM do COMSOL vêm do modelo GMD caso 3.1 (`mtl_3`, a variante mais
completa das três) — apenas o retorno à terra é, de fato, medido no COMSOL.

#### 5. Parâmetros PUL por cenário (linhas 205–229)

Para cada um dos 3 cenários (`'1'`, `'2'`, `'3'`), via `InternalPerUnitParameters` (forma
`'hybrid'`) e `PerUnitParameters`:

1. Matrizes internas núcleo+blindagem — dependem do modelo de duto, não do solo.
2. `earth_return_parameters(zg_form, yg_form)` + `quasi_tem_approx_matrices(...)` —
   **nota importante:** a formulação de `Zg`/`Yg` não modela o duto de HDPE em nenhum
   dos três cenários; a diferença entre eles está apenas no raio externo efetivo usado
   no termo próprio de imagem (que ERS/GMD estendem até `r7`). Não existe aqui uma
   "solução analítica com duto" — os três são aproximações, e a validação de fato do
   efeito do duto no retorno à terra depende da comparação com COMSOL/MATLAB pendente.

#### 6. Tabelas e validação impressas no console (linhas 231–235)

- `compare_capacitance(pul_data, c_comsol)` — tabela `C₁₁`/`C₂₂` dos três modelos de
  duto (+ COMSOL, se `cmsl_shunt_params.txt` presente).
- `validate_against_case1_reference(pul_data)` — erro relativo máximo Zs/Ysh do cenário
  `'1'` vs. MATLAB do `andreata_case1`.

#### 7. Geração de gráficos (linhas 238–245)

`SCCPlotter.compare_complete_matrices` gera, para o `key_list` atualmente ativo no
script:

| Chave (`plot_config.py`) | Conteúdo |
|---|---|
| `core_self_impedance` | `Z_cc` — autoimpedância do núcleo, fase A |
| `mutual_impedance_core_sheath` | `Z_cs` — mútua núcleo-blindagem, fase A |
| `sheath_self_impedance` | `Z_ss` — autoimpedância da blindagem, fase A |
| `earth_return_impedance_phase_a` | `Zg` — autoimpedância de retorno à terra, fase A |

Cada gráfico sobrepõe as três curvas analíticas (`DUCT_MODEL_TEMPLATE`: Underground/ERS/GMD)
e, quando disponíveis, os marcadores COMSOL/MATLAB próprios do duto
(`COMSOL_TEMPLATE`/`MATLAB_TEMPLATE`) e a referência cruzada do `andreata_case1`
(`MATLAB_CASE1_NO_DUCT_TEMPLATE`, marcador `+` laranja).

`plot_config.py` também define outras 5 chaves (`core_self_admittance`,
`mutual_admittance_core_sheath`, `sheath_self_admittance`,
`earth_return_admittance_phase_a`, `earth_return_potential_coeff_phase_a`), não incluídas
no `key_list` atual de `main()` — podem ser adicionadas à chamada
`compare_complete_matrices(key_list=[...])` conforme necessário.

`GroundReturnMTLRepresentation(__file__, mtl_0, units='centimeter').system_schematic()`
gera `Results/system_schematic.png` a partir da geometria física real (`model_0`, com
duto).

---

## Dependências

| Módulo | Papel |
|---|---|
| `models.single_core_cable.SingleCoreCableModelGenerator` | Constrói os modelos MTL (`flat_hdpe_enclosed_model`, `underground_flat_model`) a partir do JSON |
| `mtl_main.source.MulticonductorTransmissionLine` | Representa o sistema MTL multi-condutor |
| `analytical_forms.single_core_cable.InternalPerUnitParameters` | Matrizes internas do cabo (núcleo + bainha) |
| `analytical_forms.single_core_cable.PerUnitParameters` | Matrizes PUL completas com retorno à terra |
| `analytical_forms.single_core_cable.EquivalentRadiiSystems` | Raios/permissividades equivalentes ERS (área) e GMD para o duto |
| `analytical_forms.single_core_cable.apply_semiconducting_layer_correction` | Correção de camada semicondutora (núcleo e blindagem) |
| `utils.comsol_data.ComsolPostProcessor` | Leitura e processamento dos dados COMSOL (família de retorno à terra) |
| `utils.matlab_data.MatlabDataReader` | Leitura das matrizes de referência MATLAB (próprias e do `andreata_case1`) |
| `plotter.scc_plotter.SCCPlotter` | Geração dos gráficos (`compare_complete_matrices`) |
| `mtl_main.graphics.GroundReturnMTLRepresentation` | Esquema do circuito MTL |
| `plot_config.PLOT_CONFIG` | Layout, séries e limites de todos os gráficos |

---

## Limitações e Dados Pendentes

- **COMSOL do duto ainda não existe.** `get_general_parameters('cmsl_ground_return_impedance')`
  retorna `None` — bloco ignorado graciosamente. Quando adicionado, deve seguir o mesmo
  formato de colunas do `andreata_case1` (`cmsl_ground_return_impedance.txt`, sufixos
  `_rho_g_100_epsr1_1_mf_vcoil_1/2/3`); a chave de cenário `'rho_g_100_epsr1_1_mf'` já
  está fixada no código combinando com o solo do JSON (ρ=100 Ω·m, εr=1).
- **MATLAB do duto já disponível.** `MatlabDataReader(__file__, ...)` lê `.mat` em
  `Results/` com prefixo `'andreata_case2'` (ex.: `andreata_case2_frequency_range.mat`,
  `andreata_case2_series_impedance_matrix.mat` — mesmo conjunto de variáveis que o
  `andreata_case1` usa, cada caso com seu próprio prefixo `andreata_case{N}_*` para
  evitar colisão de nomes entre os diretórios `Results/`).
- **`Results/cmsl_shunt_params.txt` é do estudo antigo de cabo único** — reaproveitado
  apenas para a linha "COMSOL (energy method)" da tabela de capacitância; não é
  referência do sistema trifásico e deve ser substituído quando os dados novos chegarem
  (mesmo nome de arquivo esperado).
- **`Results/cmsl_coaxial_cable_impedance.txt`, `cmsl_series_impedance_core_excitation.txt`,
  `cmsl_series_impedance_sheath_excitation.txt`** são resíduos da linhagem antiga
  (API de cabo isolado) e não são mais lidos pelo script atual — dado morto em
  `Results/`, sem bloquear a execução.
- **Retorno à terra nunca modela o duto**, em nenhum dos 3 cenários — a única forma de
  validar o efeito real do duto no retorno à terra é a comparação COMSOL/MATLAB acima,
  ainda pendente.

Detalhes cronológicos completos das correções que levaram ao estado atual do script
(generalização para 3 fases, bug de duto concêntrico vs. excêntrico, bug de aliasing de
`center_point`, etc.) estão em `BUGS_AND_FIXES.md`.

---

## Execução

```bash
python -m testData.andreata_case2.andreata_case2
```

**Saída esperada:** tabela de capacitância e validação cruzada impressas no console, 4
gráficos (conforme `key_list` atual) mais `system_schematic.png` salvos/exibidos a
partir de `Results/`. Os dados COMSOL/MATLAB próprios do duto são opcionais — a
simulação analítica e a validação cruzada com o `andreata_case1` rodam normalmente sem
eles.

---

## Referências

- **Andreata, L. E. B.** (2025). *Análise das Características de Propagação e de
  Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não Metálicos
  no Contexto de Parques Eólicos.* UFMG. https://hdl.handle.net/1843/2086
- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return
  Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de
  Montréal. (formulações `magalhaes_xue` de retorno à terra)
- **Lafaia, I.** (2015). Método GMD para permissividade equivalente de dutos (caso 3.1).
