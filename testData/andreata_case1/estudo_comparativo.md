# Estudo comparativo — andreata_case1 / case2 / case3

## 1. Arquitetura comum (o "esqueleto" dos três `andreata_caseN.py`)

Os três scripts seguem exatamente o mesmo molde de pipeline:

```
imports (case_utils, ComsolPostProcessor, MatlabDataReader, SCCPlotter,
         SingleCoreCableModelGenerator, MulticonductorTransmissionLine,
         InternalPerUnitParameters/PerUnitParameters) + PLOT_CONFIG
  ↓
cable_generator = SingleCoreCableModelGenerator(__file__)
model = cable_generator.<método específico do caso>()
model = apply_semiconducting_layer_correction(model, core, sheath)
mtl_model_a = MulticonductorTransmissionLine(model)
  ↓
pul_data = {'frequencies': logspace(-2,7,90), 'comsol': {...}, 'scenarios': {...}}
  ↓
[opcional] ComsolPostProcessor(__file__) → earth_return + quasi_tem por cenário
MatlabDataReader(__file__).get_scc_scenario_data(prefix, conductor_order)
  ↓
InternalPerUnitParameters + PerUnitParameters por cenário (loop)
  ↓
SCCPlotter(__file__, pul_data, PLOT_CONFIG).compare_internal_matrices(...)
GroundReturnMTLRepresentation(...).system_schematic(); plt.show()
```

Essa uniformidade é resultado direto do trabalho de reescrita já registrado nos `BUGS_AND_FIXES.md` — `andreata_case2.py` foi deliberadamente reescrito "no molde do case1/case3" (Item 5 do seu changelog), e é isso que se confirma lendo o código hoje.

## 2. Tabela comparativa

| Aspecto | case1 | case2 | case3 |
|---|---|---|---|
| Geometria | 3 SCC enterrados direto (`underground_flat_model`) | 3 SCC, cada um em duto HDPE individual (`flat_hdpe_enclosed_model`) | 3 SCC + 1 ECC, sem duto (`flat_scc_with_ecc_cable_model`) |
| Nº condutores | 6 (homogêneo) | 6 (homogêneo) | 7 (heterogêneo: 3×2 + 1×1) |
| Eixo de comparação | formulação de solo (`magalhaes_xue` × `deconti`) × 3 perfis (ρ, εr) = 6 cenários | modelo de duto (Underground / ERS / GMD) = 3 cenários, 1 só formulação de solo | mesmo padrão do case1: 6 cenários (formulação × perfil de solo) |
| `ComsolPostProcessor` chamado? | **sim** (retorno à terra + internal combinado) | **sim** (idêntico ao case1) | **não** — bloco inteiro ausente |
| `MatlabDataReader` | 1 instância, `conductor_order=[0,3,1,4,2,5]` | **2 instâncias** — a própria + uma apontando para `Results` do case1 (validação cruzada) | 1 instância, `conductor_order=[0,3,1,4,2,5,6]` |
| Funções auxiliares locais | nenhuma | `_override_sheath_insulation`, `compare_capacitance`, `validate_against_case1_reference` | nenhuma |
| `plot_config.py` | 328 linhas | 333 linhas | 352 linhas |

## 3. Redundâncias verificadas (com localização exata)

**3.1 — `plot_config.py`: 131 linhas idênticas byte-a-byte entre case1 e case3**
Confirmado por `diff`: as linhas 1–131 de `testData/andreata_case1/plot_config.py` e `testData/andreata_case3/plot_config.py` são **exatamente iguais** — `COMSOL_TEMPLATE`, `XUE_TEMPLATE`, `VANCE_TEMPLATE`, `MATLAB_TEMPLATE`, `DECONTI_TEMPLATE` e os seis `PLOT_TPL_*`. O case2 duplica um subconjunto disso (`MATLAB_TEMPLATE` + os 6 `PLOT_TPL_*`, ~53 linhas) verbatim.

**3.2 — `VANCE_TEMPLATE` é código morto em dois arquivos**
Definido em `andreata_case1/plot_config.py:40` e `andreata_case3/plot_config.py:40`, mas nunca referenciado em nenhum `PLOT_CONFIG` dos três casos (confirmado via grep). Sobra de copy-paste de outro caso que usa formulação de Vance.

**3.3 — `utils/case_utils.py`: função duplicada + função quebrada, ambas não usadas pelos três casos**
- `matrix_to_string` está definida **duas vezes** (`case_utils.py:35` e `case_utils.py:108`) — a segunda sobrescreve a primeira silenciosamente (dead code, ~20 linhas).
- `print_real_matrix` (`case_utils.py:116`) chama `_matrix_to_string(...)` — nome que **não existe** no módulo (a função real chama-se `matrix_to_string`, sem underscore). Se algo chegar a chamar `print_real_matrix`, estoura `NameError`. Como `from utils.case_utils import *` é usado nos três `andreata_caseN.py`, o erro entraria silenciosamente no namespace de cada um sem ser notado — só dispararia se alguém tentasse usar a função.
- Confirmado por grep: nenhum dos três casos `andreata_*` usa `print_real_matrix`, `verify_kelvin_functions`, `matrix_viewer`, `format_complex_number` ou `format_scientific_notation` — são utilitários de outras linhagens de teste (`mom_models.py`, `deConti_models.py`, casos `_s40`/`_s50`/`_s100`) importados de graça via `import *`.

**3.4 — Blocos de comentário quase idênticos (copy-paste) em `andreata_case1.py` e `andreata_case2.py`**
O bloco "COMSOL de impedância interna (medição legada...)" em `andreata_case1.py:106-114` e o equivalente em `andreata_case2.py:205-216` são o mesmo texto com pequenas adaptações — inclusive repetindo a mesma nota de risco ("compartilha `pul_data['comsol']['frequencies']`... hoje inofensivo... se os dois vierem a coexistir"). Sinaliza que essa lógica (carregar `cmsl_internal_impedance_matrix.txt` combinado) é candidata natural a virar um método único em `ComsolPostProcessor` em vez de bloco replicado no script de cada caso.

## 4. Divergência estrutural mais relevante: case3 não carrega COMSOL

`andreata_case3.py` **não instancia `ComsolPostProcessor` em nenhum momento** — o bloco inteiro que case1/case2 têm (`"Construindo matrizes COMSOL..."` + `"Carregando dados COMSOL de impedância interna..."`, ~35 linhas) simplesmente não existe no case3. Isso é coerente com o fato de não haver `Results/cmsl_ground_return_impedance.txt` para esse caso (documentado no Item 9, "pendente", de `BUGS_AND_FIXES.md` do case3).

O problema é que **`andreata_case3/plot_config.py`** já vem com `comsol_series_to_plot: COMSOL_TEMPLATE` configurado em 5 gráficos (`mutual_impedance_phase_a_sheath_ecc`, `self_impedance_ecc`, `self_admittance_ecc`, `earth_return_impedance_ecc`, `earth_propagation_constant`) — configuração morta hoje (nunca populada, `self.cmsl.get('frequencies')` sempre `None` → o guard em `scc_plotter.py` pula silenciosamente), mas que vai **quebrar exatamente como os Itens 4/5 já quebraram** no dia em que alguém adicionar o arquivo COMSOL: o Item 9 já documenta que `get_quasi_tem_approx_matrices` (`utils/comsol_data.py:731`, `M = 2` hardcoded na linha 738) não tem a generalização por `block_sizes` que o lado analítico ganhou no Item 5 — vai quebrar em `np.kron(4×4, ones((2,2)))` = 8×8 contra uma `Zi` 7×7.

## 5. Oportunidades de limpeza/otimização, priorizadas

| # | Ação | Ganho |
|---|---|---|
| 1 | Extrair `COMSOL_TEMPLATE`/`XUE_TEMPLATE`/`VANCE_TEMPLATE`/`MATLAB_TEMPLATE`/`DECONTI_TEMPLATE`/`PLOT_TPL_*` para um módulo compartilhado (`testData/andreata_common/plot_templates.py` ou similar) importado pelos três `plot_config.py` | Elimina ~130 linhas duplicadas × 2 arquivos; qualquer ajuste de estilo (cor, marcador, escala) passa a valer para os três de uma vez, sem risco de um caso ficar dessincronizado |
| 2 | Remover `VANCE_TEMPLATE` dos dois arquivos onde está morto | Reduz ruído; se for reintroduzido depois, adicionar junto do uso real |
| 3 | Corrigir/remover `case_utils.py:108` (segunda definição de `matrix_to_string`) e `print_real_matrix` (referência quebrada a `_matrix_to_string`) | Elimina dead code + um bug latente que hoje só não estourou porque nada chama a função |
| 4 | Resolver o Item 9 pendente: passar `block_sizes` explícito para `ComsolPostProcessor.get_quasi_tem_approx_matrices`, reusando `_expand_by_block_sizes` (já existe, criada para o Item 5) | Fecha a lacuna antes que alguém adicione dado COMSOL a case3 e recaia no mesmo bug já corrigido no lado MATLAB/analítico |
| 5 | Extrair o bloco "COMSOL de impedância interna combinada" (idêntico em case1.py/case2.py) para um método de conveniência em `ComsolPostProcessor`, ex. `load_internal_and_earth_return(cmsl_processor, mtl_model, scenarios)` | Reduz ~35 linhas duplicadas por script para uma chamada; comentário de risco fica documentado uma vez só, no código, não em dois lugares |
| 6 | Se/quando case3 ganhar dado COMSOL de retorno à terra, copiar o bloco de carregamento de case1.py — hoje o `plot_config.py` já está pronto, só falta o loader | Evita que o "pendente" vire outro ciclo de debug como os Itens 4/5/6/7/8 |

## 6. O que **não** deveria ser mexido sem cuidado

- As diferenças reais de arquitetura (case2 com `_override_sheath_insulation`/`compare_capacitance`/`validate_against_case1_reference`, a segunda instância de `MatlabDataReader` apontando pro case1) refletem necessidade física genuína — case2 é o único com modelo de duto e validação cruzada contra case1 sem duto — não são redundância a eliminar, mas também não valem a pena generalizar para os outros dois casos sem necessidade concreta.
- `MergedComsolDataReader` (em `utils/comsol_data.py`) parece redundante com `ComsolDataReader._parse_single_file`, mas na verdade serve uma linhagem diferente (`testData/bare_wire`) — não é usada por nenhum dos três casos andreata, então não é redundância *dentro* do escopo pedido, só duas estratégias de parsing coexistindo no mesmo módulo para casos diferentes.

## 7. Plano de otimização — execução e resultados

As seis ações da seção 5 foram implementadas em seis fases sequenciais, cada uma testada (import + execução ponta-a-ponta dos três `andreata_caseN.py`, mais regressão em casos vizinhos) antes de avançar para a próxima. Duas decisões de projeto foram definidas antes de começar:

- **Local do módulo de templates compartilhado:** `testData/andreata_common/plot_templates.py` — pasta nova e neutra, em vez de `utils/` (que mistura infraestrutura de dados com config de apresentação) ou de dentro de um dos casos (que criaria uma dependência estranha de um caso "irmão").
- **`print_real_matrix` quebrada:** corrigida (não removida) — mantém a função disponível para uso futuro.

### Fase 0 — Limpeza de dead code

- `utils/case_utils.py`: removida a segunda definição duplicada de `matrix_to_string` (linha 108, sobrescrevia a primeira); corrigida a chamada quebrada em `print_real_matrix` (`_matrix_to_string` → `matrix_to_string`).
- `andreata_case1/plot_config.py` e `andreata_case3/plot_config.py`: removido `VANCE_TEMPLATE` (código morto).

**Validação:** `print_real_matrix` testada isoladamente (antes estourava `NameError`); os três `andreata_caseN.py` rodados de ponta a ponta sem erros novos.

### Fase 1 — Templates de `plot_config.py` extraídos para módulo compartilhado

- Criado `testData/andreata_common/plot_templates.py` com duas camadas: base comum aos três (`MATLAB_TEMPLATE`, `PLOT_TPL_REAL/IMAG/RESISTANCE/INDUCTANCE/CONDUCTANCE/CAPACITANCE`) e camada específica de case1/case3 (`COMSOL_TEMPLATE` de 3 perfis, `XUE_TEMPLATE`, `DECONTI_TEMPLATE`).
- `andreata_case1/plot_config.py` e `andreata_case3/plot_config.py`: as 131 linhas duplicadas viraram um import de 5 linhas.
- `andreata_case2/plot_config.py`: importa só a camada base; mantém local o que é específico dele (`COMSOL_TEMPLATE` de 1 perfil, `MATLAB_CASE1_NO_DUCT_TEMPLATE`, `DUCT_MODEL_TEMPLATE`, `INTERNAL_COMSOL_TEMPLATE`).

**Validação:** checagem estrutural (`PLOT_CONFIG` dos três casos com o mesmo número de chaves de antes, templates importados resolvidos por identidade/igualdade dentro dos dicts) + execução ponta-a-ponta dos três casos, mesma saída de console de antes.

### Fase 2 — Bloco de carregamento COMSOL duplicado extraído

- Novo método `ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios()` (`utils/comsol_data.py`), encapsulando o bloco antes duplicado em `andreata_case1.py:87-123` e `andreata_case2.py:183-224` (retorno à terra por cenário + impedância interna combinada). O comentário de risco sobre `pul_data['comsol']['frequencies']` compartilhado passou a viver uma vez só, no docstring do método.
- `andreata_case1.py` e `andreata_case2.py` reescritos para chamar o novo método em vez do bloco inline.

**Validação:** os três `andreata_caseN.py` + `scc_flat_xue.py` (outro consumidor de `ComsolPostProcessor`, não tocado) rodados de ponta a ponta; mesmos números de saída de case2 (tabela de capacitância, erros de validação cruzada 0.06%/0.15%). Achado à parte, fora do escopo: `testData/scc_flat_xue/scc_flat_xue.py:129` chama `plotter.scc_earth_return_impedance_matrix()`, método inexistente em `SCCPlotter` — bug pré-existente, não introduzido por este trabalho, não corrigido.

### Fase 3 — Item 9 (pendente) resolvido: `get_quasi_tem_approx_matrices` generalizado por `block_sizes`

- `utils/comsol_data.py`: `M = 2` hardcoded + laço `np.kron` por frequência trocados por `block_sizes = internal_matrices['block_sizes']` + `_expand_by_block_sizes` (reaproveitada de `analytical_forms/single_core_cable.py`, mesma função criada para o Item 5 do lado analítico).

**Validação:**
- Com dados reais (`scc_flat_xue`, caso homogêneo `block_sizes=[2,2,2]`): resultado **bit-idêntico** (`np.array_equal`) ao `np.kron` antigo, para `Zg` e `Pg`.
- Sinteticamente (`block_sizes=[2,2,2,1]`, formato do caso ECC): monta corretamente uma matriz 7×7 com os blocos no lugar certo — o código antigo quebraria com `np.kron(4×4, ones((2,2)))` → 8×8, incompatível com a `Zi` 7×7.
- Regressão nos três casos + `scc_flat_xue`: comportamento idêntico ao anterior.

### Fase 4 — COMSOL habilitado (guardado) em `andreata_case3.py`

- Adicionado o bloco `ComsolPostProcessor(__file__)` + `load_scc_earth_return_and_internal_scenarios(...)` a `andreata_case3.py`, ausente até então. Como não existe `Results/cmsl_ground_return_impedance.txt` nem `cmsl_internal_impedance_matrix.txt` para este caso, o bloco fica inerte hoje (só emite os dois avisos de dado ausente), mas já protegido pelo fix da Fase 3 e pronto para os 5 gráficos que `plot_config.py` já tem configurados com `comsol_series_to_plot`.

**Validação:** `andreata_case3.py` rodado de ponta a ponta — exit 0, dois avisos esperados, os mesmos 6 cenários calculados, sem warnings novos.

### Fase 5 — Regressão ampla

Rodados `scc_34kV_andreata`, `hdpe_2000mm2`, `hdpe_300mm2`, `hdpe_ecc_2000mm2`, `scc_138kV_prysmian` (além dos três andreata e do `scc_flat_xue` já cobertos nas fases anteriores) — todos com exit 0, sem tracebacks ou erros ocultos no log.

### Estado final

Todas as seis ações da seção 5 foram concluídas e validadas. Nenhuma mudança de comportamento numérico foi introduzida nos casos homogêneos (case1, case2, `scc_flat_xue`); o caso heterogêneo (case3/ECC) ganhou proteção contra um bug que ainda não havia sido exercitado, mas que quebraria assim que dado COMSOL de retorno à terra fosse adicionado a esse caso.
