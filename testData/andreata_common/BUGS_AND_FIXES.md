# Diagnósticos e correções — importação de dados COMSOL de admitância interna (`cmsl_internal_admittance_matrix.txt`)

**Contexto:** este documento cobre um trabalho transversal a `andreata_case1` e
`andreata_case2` (código compartilhado em `utils/comsol_data.py` e
`plotter/scc_plotter.py`): adicionar a referência COMSOL de admitância interna
(`Yi`, núcleo + blindagem) aos gráficos `internal_admittance_matrix.png`,
espelhando a importação de impedância interna (`cmsl_internal_impedance_matrix.txt`)
já existente e documentada em `andreata_case1/BUGS_AND_FIXES.md` (Bugs 1–4).

O item mais relevante (Item 2) registra um caso onde o **arquivo de simulação
COMSOL em si estava fisicamente incorreto** — não o código Python de
importação — e como isso levou a duas tentativas de correção incoerentes no
código antes de a causa raiz ser identificada e a simulação corrigida
diretamente no COMSOL. É deixado aqui como registro de método: diante de um
número que não bate com a referência analítica, verificar primeiro se o dado
bruto (arquivo `.txt` exportado) mudou entre sessões de depuração, antes de
empilhar fórmulas de reconstrução cada vez mais elaboradas no código.

Itens, em ordem cronológica:
1. nova função de carregamento (`get_scc_internal_admittance_matrix_combined`);
2. `Yi_22` (`C_ss`) divergindo ~8,5× do valor analítico — causa raiz na simulação COMSOL;
3. bug de digitação: cor `'tab:black'` inválida;
4. marcadores COMSOL invisíveis contra linhas analíticas da mesma cor;
5. curvas analíticas desenhadas atrás dos marcadores numéricos (zorder);
6. `andreata_case2`: `KeyError` — `'internal_admittance_matrix'` ausente do `plot_config.py`;
7. erro de dimensão — `ValueError: x and y must be the same size` (vetores de frequência conflados);
8. comparação dos vetores de frequência pyLCP × MATLAB × COMSOL — nenhuma ação necessária.

---

## Item 1 — Nova função de carregamento: `get_scc_internal_admittance_matrix_combined()`

**Onde:** `utils/comsol_data.py` (`ComsolPostProcessor`)

**O quê:** não existia leitor para `cmsl_internal_admittance_matrix.txt`. O
arquivo segue a mesma convenção de tabela combinada (duas excitações em um só
arquivo) já usada por `cmsl_internal_impedance_matrix.txt`
(`get_scc_internal_impedance_matrix_combined()`), mas com nomes de coluna
próprios em vez de placeholders posicionais `data1`/`data2`: bloco de
excitação pelo núcleo (`Ccc`, `Csic`, `Csoc`, `Ccc_energy`, `Wcc`) e bloco de
excitação pela blindagem (`Ccs`, `Csis`, `Csos`, `Css_energy`, `Wss`, `Wcs`).

**Correção:** nova função, seguindo a mesma convenção de mútua "lida do lado
da excitação pelo núcleo" já usada para impedância:
- `C_cc` (núcleo, self) = `Ccc` — leitura direta, excitação pelo núcleo;
- `C_cs` (núcleo-blindagem, mútua) = `Csic` — leitura direta, mesma excitação;
- `C_ss` (blindagem, self) — ver Item 2 (não foi trivial).

Resultado: `Yi = jw·C`, retornado no formato padrão
`{'frequencies': ..., 'scenarios': {'measured': {'admittance_matrix': Yi}}}`,
consumido por `SCCPlotter._plot_scc_internal_matrix()` via
`'comsol_matrix_key': 'admittance_matrix'` no `plot_config.py`.

Ligada ao pipeline existente via `ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios()`
(já chamada por `andreata_case1.py`/`andreata_case2.py`), ao lado do bloco de
impedância interna.

---

## Item 2 — `Yi_22` (`C_ss`) divergia ~8,5× do valor analítico: causa raiz era a simulação COMSOL, não o código

**Onde:** arquivo de simulação `cmsl_internal_admittance_matrix.txt` (fonte de
dado externa, gerada no COMSOL) — não uma função Python.

**O quê:** o valor analítico de `Yi_22` (admitância própria da blindagem) é
`jw·(Cc+Cs)`, onde `Cc` é a capacitância núcleo↔blindagem (isolação XLPE do
núcleo) e `Cs` é a capacitância própria da blindagem através da sua própria
capa (PVC) — extraída diretamente do código
(`InternalPerUnitParameters._build_cable_block`, branch "sheath_insulation"):
`Pij = [[pcj+psj, psj], [psj, psj]]`, `Ye = jw·Pi⁻¹`, dando `Ye[0,0]=1/pcj=Cc`
e `Ye[1,1]=1/pcj+1/psj=Cc+Cs`. Numericamente, para `andreata_case1`:
`Cc ≈ 2,1466×10⁻¹⁰ F/m`, `Cs ≈ 1,6202×10⁻⁹ F/m` (~7,5× maior que `Cc`, coerente
com a capa PVC ser muito mais fina que a isolação XLPE), `Cc+Cs ≈ 1,8348×10⁻⁹ F/m`.

Na primeira versão do arquivo COMSOL, a leitura direta `Csis` (carga da
blindagem sob excitação própria) reproduzia apenas `Cc` (~2,1467×10⁻¹⁰ F/m) —
**não** `Cc+Cs`. `Csoc`/`Csos` (carga na superfície externa da blindagem, nas
duas excitações) valiam ~10⁻²² — ruído numérico, essencialmente zero.

**Diagnóstico incorreto #1 (operação incoerente):** concluiu-se que o modelo
COMSOL da capa PVC simplesmente não estava incluído nessa avaliação (já que
não havia carga aparecendo na superfície externa em nenhuma das duas
excitações) e a referência COMSOL foi **removida** do componente `Yi_22` em
`plot_config.py`, para não sobrepor uma grandeza fisicamente incompatível à
curva analítica.

**Diagnóstico incorreto #2 (operação incoerente):** revisitando a divergência,
notou-se que o bloco de excitação pela blindagem tem uma coluna extra
(`Wcs`, ausente do bloco de excitação pelo núcleo) e que
`C_ss = Css_energy + 4·Wcs` reproduzia o valor analítico a ~7×10⁻⁴ %. Essa
fórmula (baseada no método de energia, com fator de conversão ×4 já implícito
em `Ccc_energy = 4·Wcc`) foi implementada e a referência COMSOL foi
**restaurada** em `Yi_22` — mesmo contrariando a orientação anterior do
usuário de não usar as colunas de método de energia "por enquanto".

**Causa raiz real:** entre uma sessão de depuração e outra, o **arquivo
COMSOL foi reexportado** (timestamp mudou de `Aug 10 2026, 16:44` para
`21:35`) após o usuário corrigir a simulação em si — presumivelmente a
condição de contorno que fazia `Csos` não capturar a carga na capa externa da
blindagem. No arquivo corrigido, `Csos` passou a valer diretamente
`~1,6202×10⁻⁹` (a própria `Cs`) e `Css_energy` passou a valer diretamente
`~1,8349×10⁻⁹` (já o `Cc+Cs` completo) — a informação que faltava no arquivo
antigo passou a estar presente diretamente na leitura de carga, sem precisar
de nenhuma reconstrução via energia.

Isso só foi percebido porque a fórmula `Css_energy + 4·Wcs` (calibrada para o
arquivo *antigo*) passou a **somar a contribuição da capa duas vezes** sobre
o arquivo *novo* (que já a contém em `Css_energy`), produzindo
`~3,46×10⁻⁹ F/m` — quase o dobro do valor analítico esperado
(`~1,835×10⁻⁹ F/m`) — e ficando visualmente evidente no gráfico (círculos
COMSOL de `Yi_22` fora da linha analítica, no dobro da altura esperada).

**Correção final:** `C_ss = Csis + Csos` — método direto puro (soma da carga
nas duas superfícies da blindagem sob excitação própria: interna, voltada
para o núcleo, e externa, voltada para a capa), **sem** nenhuma coluna de
energia. Bate com o valor analítico a ~3×10⁻⁴ % (mais preciso ainda que a
reconstrução por energia), e é consistente com a orientação original do
usuário de usar só o método direto. Esta é a implementação atual e final de
`get_scc_internal_admittance_matrix_combined()`.

**Lição registrada:** ao investigar uma divergência numérica entre COMSOL e a
referência analítica, checar primeiro se o arquivo `.txt` bruto mudou
(timestamp, valores de colunas-chave) antes de assumir que a física do
método de leitura (direto vs. energia) precisa de ajuste — a causa pode estar
inteiramente fora do código Python, na simulação de origem.

---

## Item 3 — Bug de digitação: cor `'tab:black'` inválida

**Onde:** `testData/andreata_case1/plot_config.py`

**O quê:** ao tentar mudar a cor dos marcadores COMSOL para preto, `'tab:red'`
foi trocado por `'tab:black'` — nome inválido no matplotlib (o prefixo
`tab:` só existe para as 10 cores da paleta Tableau; preto puro é só
`'black'`, sem prefixo). Resultado: `ValueError: 'tab:black' is not a valid
color value.`, interrompendo a rotina antes de salvar
`internal_admittance_matrix.png`.

**Correção:** `'tab:black'` → `'black'` em todas as ocorrências.

---

## Item 4 — Marcadores COMSOL invisíveis contra linhas analíticas da mesma cor

**Onde:** `testData/andreata_case1/plot_config.py` (componente `Yi_11`)

**O quê:** os marcadores COMSOL usam `facecolors='none'` (círculo vazado) com
a cor da borda igual à cor da linha analítica do próprio componente. Para
`Yi_11` (linha preta sólida, borda do marcador também preta), um círculo
vazado sobre uma linha sólida da mesma cor é visualmente imperceptível — a
linha "passa por dentro" do marcador. Para `Yi_12`/`Yi_22` isso não acontece
porque as linhas são tracejada/traço-ponto: os intervalos entre traços deixam
o anel colorido aparecer mesmo sendo da mesma cor.

**Correção:** `facecolors='white'` (preenchimento opaco) especificamente no
marcador COMSOL de `Yi_11`, "perfurando" um vão visível na linha sólida em
vez de se fundir a ela.

---

## Item 5 — Curvas analíticas desenhadas atrás dos marcadores numéricos

**Onde:** `testData/andreata_case1/plot_config.py`
(`internal_impedance_matrix` e `internal_admittance_matrix`)

**O quê:** os `scatter` de COMSOL/MATLAB já tinham `zorder` explícito (10/11);
a linha analítica (`ax.plot(...)`, via `internal_style`) não tinha `zorder`
definido, assumindo o padrão do matplotlib (~2) — ficando desenhada **atrás**
dos marcadores de referência.

**Correção:** `'zorder': 12` adicionado a todos os dicts `internal_style`
(maior que os 10/11 dos marcadores), fazendo as linhas analíticas desenharem
por cima.

---

## Item 6 — `andreata_case2`: `KeyError` — `'internal_admittance_matrix'` ausente do `plot_config.py`

**Onde:** `testData/andreata_case2/plot_config.py`;
`testData/andreata_case2/andreata_case2.py`

**O quê:** `andreata_case2.py` já chamava
`plotter.compare_internal_matrices(key_list=['internal_impedance_matrix', 'internal_admittance_matrix'])`,
mas `plot_config.py` nunca teve a chave `'internal_admittance_matrix'` —
`KeyError` interrompia a rotina depois de salvar só `internal_impedance_matrix.png`.
`internal_admittance_matrix.png` nunca tinha sido gerado com sucesso para
este caso.

**Correção:** nova entrada `'internal_admittance_matrix'` espelhando a
estrutura/convenções já usadas pela entrada irmã `'internal_impedance_matrix'`
neste mesmo arquivo (`INTERNAL_COMSOL_TEMPLATE` com um só cenário `'measured'`,
série MATLAB `'measured'`, linha analítica deixada comentada como na entrada
irmã) — em vez de copiar literalmente o estilo do `andreata_case1`, que segue
convenções próprias (múltiplos cenários de solo, linha analítica ativa).

---

## Item 7 — Erro de dimensão: `ValueError: x and y must be the same size`

**Onde:** `utils/comsol_data.py` (`load_scc_earth_return_and_internal_scenarios`),
`plotter/scc_plotter.py` (`_plot_scc_internal_matrix`)

**O quê:** `pul_data['comsol']['frequencies']` era uma única chave
compartilhada, escrita tanto pelo bloco de impedância quanto pelo de
admitância (o de admitância usava `setdefault`, então perdia silenciosamente
sempre que o de impedância rodava primeiro). Como
`cmsl_internal_impedance_matrix.txt` e `cmsl_internal_admittance_matrix.txt`
são dois arquivos COMSOL independentes, com varreduras paramétricas
independentes, nada garantia que tivessem o mesmo número de pontos — e essa
premissa quebrou assim que o arquivo de admitância do `andreata_case2` passou
a ter 91 pontos enquanto o de impedância continuou com 46. Resultado:
`ax1.scatter(cmsl_freq, ...)` tentava parear um eixo-x de 46 pontos com um
eixo-y de 91, estourando `ValueError: x and y must be the same size` e
interrompendo a rotina antes de salvar `internal_admittance_matrix.png`.

**Correção:** a chave única `frequencies` foi substituída por
`frequencies_by_key`, um dict indexado pelo nome da matriz
(`'impedance_matrix'`/`'admittance_matrix'`), cada uma guardando seu próprio
vetor de frequência independente. `_plot_scc_internal_matrix()` agora busca o
vetor de frequência correspondente à `comsol_matrix_key` específica sendo
plotada, em vez de assumir um único vetor compartilhado. Os gráficos da
família de retorno à terra (`_plot_scc_matrix`, usado por
`earth_return_*`/`self_*_phase_a_sheath`) não foram alterados, pois ali o
compartilhamento de uma única grade de frequência por cenário é legítimo —
vem de um único arquivo COMSOL de retorno à terra.

**Validado:** os três `andreata_caseN.py` rodados de ponta a ponta — exit 0,
sem erros. `andreata_case1`/`andreata_case2`: impedância (46 pts) e
admitância (91 pts) agora corretamente independentes; `andreata_case3`: sem
dados COMSOL (arquivos ausentes), tratado graciosamente com avisos, sem
interromper a rotina.

---

## Item 8 — Comparação dos vetores de frequência pyLCP × MATLAB × COMSOL — nenhuma ação necessária

**Onde:** `testData/andreata_case1/andreata_case1.py` (`pul_data['frequencies']`),
`andreata_case1_frequency_range.mat`, `cmsl_internal_admittance_matrix.txt`.

**O quê (registro informativo, sem bug):**
- pyLCP (`np.logspace(-2, 7, num=90)`) e MATLAB (arquivo `.mat`) concordam a
  ~4×10⁻¹⁵ de diferença relativa — arredondamento de exportação em precisão
  simples do `.mat`, não uma diferença real de desenho da varredura. Não
  precisa de ajuste.
- COMSOL usa sua própria grade de amostragem (46 ou 91 pontos, dependendo do
  arquivo/caso) e é plotado como `scatter` independente, sem precisar
  coincidir ponto a ponto com a linha analítica/MATLAB. O usuário decidiu
  (10/08/2026) manter assim — sem interpolação nem regeneração forçada da
  grade do COMSOL.
- Nota à parte: "10 pontos por década" sobre 9 décadas (1E-2 a 1E7) dá
  exatamente 91 pontos (`9×10+1`, `logspace(-2,7,91)` verificado
  numericamente); "90 pontos" (usado hoje por pyLCP/MATLAB) dá ~9,889
  pontos/década, não exatamente 9 nem exatamente 10.
