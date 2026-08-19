# Diagnósticos e correções — `andreata_case4` (Configuração 4, Figura 5.4)

**Contexto:** `andreata_case4` combina o duto HDPE do `andreata_case2` com o
ECC heterogêneo do `andreata_case3`. Sua criação exigiu um novo método
gerador (`flat_hdpe_enclosed_with_shared_ecc_model`) e um novo `mtl_type`
(`scc-flat-hdpe-ecc`), feitos junto com a padronização do campo `"type"` dos
quatro JSONs andreata (`"scc-flat"`, `"scc-flat-hdpe"`, `"scc-flat-ecc"`,
`"scc-flat-hdpe-ecc"`). Três problemas reais surgiram nesse trabalho, todos
descobertos e corrigidos ainda na primeira implementação.

---

## Item 1 — `'HDPE'` (maiúsculo) não registrado após tornar `"type"` funcional

**Onde:** `mtl_main/strategy.py` (`mtl_strategy_factory`), `mtl_main/graphics.py`.

**Contexto:** ao tornar `underground_flat_model()`/`flat_hdpe_enclosed_model()`
sensíveis ao campo `"type"` do JSON (antes hardcodavam `'scc'`/`'hdpe'`,
ignorando o JSON por completo), confirmou-se por grep que nenhum dos 8 casos
irmãos que compartilham esses dois métodos tinha campo `"type"` no nível
raiz do JSON — **exceto** `hdpe_300mm2.json`/`hdpe_2000mm2.json`, que já
tinham `"type": "HDPE"` (maiúsculo, não `"hdpe"`) havia tempo, sem que isso
importasse (o hardcode anterior o ignorava). Ao tornar o campo funcional,
esse valor passou a vazar para `model['type']`, e `'HDPE'` não estava
registrado em nenhum lugar — `mtl_strategy_factory` lançava
`ValueError: Unknown or unsupported MTL type: HDPE`.

**Causa raiz:** a verificação inicial de segurança (grep pela string
literal `"type"` nos 8 JSONs, usando um glob com chaves `{a,b,c}/*.json`)
retornou "nenhum resultado" incorretamente — o glob com chaves não expandiu
como esperado na ferramenta de busca usada, mascarando a presença de
`hdpe_300mm2.json`/`hdpe_2000mm2.json` na lista. Uma segunda verificação,
usando o parser JSON do Python diretamente arquivo a arquivo (mais lenta,
porém confiável), revelou o valor `"HDPE"` nos dois arquivos.

**Correção:** `'HDPE'` registrado como alias explícito em
`mtl_strategy_factory` (→ `SingleCoreCableInHDPEStrategy`, mesma classe de
`'hdpe'`) e nos mesmos pontos de `mtl_main/graphics.py` onde `'hdpe'` já
aparecia (branch de título/`h_factor` e tupla de exceção de
`depth_ref_conductor`).

**Validado:** `hdpe_300mm2.py`/`hdpe_2000mm2.py` voltaram a rodar até o fim
sem erro, com os mesmos gráficos/esquemáticos de antes.

**Lição:** ao verificar "nenhum outro caso usa X" via busca textual antes de
uma mudança que depende disso, preferir um parser real (ex.: `json.load`)
a um grep com glob complexo — um falso negativo aqui quase chegou a
implementação sem essa rede de segurança.

---

## Item 2 — Dispatch heterogêneo restrito à string exata `'scc-flat-ecc'`

**Onde:** `analytical_forms/single_core_cable.py:944`
(`InternalPerUnitParameters.matrices`).

**Sintoma original:** ao rodar `andreata_case4.py` pela primeira vez (após
corrigir o Item 1), o script chegava a `EquivalentRadiiSystems(mtl_ers_base)`
e quebrava com `KeyError: 'hdpe'` em
`analytical_forms/single_core_cable.py:265` — mas o erro real não estava
nessa linha, e sim numa causa anterior:

**Causa raiz (duas camadas):**
1. **`model_ers_base`** é construído via `flat_hdpe_enclosed_model()` sobre
   o `SingleCoreCableModelGenerator` de `andreata_case4.json` — cujo `"type"`
   de nível raiz é `"scc-flat-hdpe-ecc"` (descreve o caso como um todo, para
   o modelo heterogêneo). Como `flat_hdpe_enclosed_model()` agora lê
   `self.input_data.get('type', 'hdpe')` (Item 0 da padronização), e o JSON
   *tem* um `"type"`, o modelo homogêneo resultante herdou `type =
   'scc-flat-hdpe-ecc'` em vez de `'hdpe'` — mapeando para
   `SingleCoreCableWithECCStrategy` (heterogênea) em vez de
   `SingleCoreCableInHDPEStrategy`. O `context.scc` heterogêneo (agrupado
   por cabo) não tem chave `'hdpe'`, daí o `KeyError`.
   **Corrigido** forçando `model_ers_base['type'] = 'hdpe'` explicitamente
   em `andreata_case4.py`, logo após a geração — este modelo é sempre
   homogêneo por construção (nunca tem ECC), então o `"type"` do JSON (que
   descreve o caso heterogêneo como um todo) não se aplica a ele.
2. **Consequência mais séria, só percebida ao investigar a primeira:** pelo
   mesmo motivo, `model_1`/`model_2`/`model_3` (construídos via
   `flat_scc_with_ecc_cable_model()`, usada para os 3 cenários analíticos)
   *também* herdam `type = 'scc-flat-hdpe-ecc'` do JSON, em vez do
   `'scc-flat-ecc'` que esse método usa como default. Isso não quebra o
   registro de estratégia (ambos os tipos mapeiam para
   `SingleCoreCableWithECCStrategy` em `mtl_strategy_factory` — deliberado,
   ver Item 2 de "Mudanças" no plano), mas quebraria silenciosamente
   `InternalPerUnitParameters.matrices()`: o dispatch para o caminho
   heterogêneo (`_matrices_heterogeneous`, que sabe montar blocos
   `[2,2,2,1]`) checava a string exata `mtl_type == 'scc-flat-ecc'` — com
   `mtl_type == 'scc-flat-hdpe-ecc'`, cairia no ramo genérico/homogêneo,
   incompatível com o `context.scc` agrupado por cabo (mesma classe de erro
   do Item 4 de `andreata_case3/BUGS_AND_FIXES.md`, agora por um caminho
   diferente).

**Correção:** generalizado o dispatch em
`analytical_forms/single_core_cable.py:944` para
`if self.model.mtl_type in ('scc-flat-ecc', 'scc-flat-hdpe-ecc'):` — a
verificação correta é "este `mtl_type` mapeia para
`SingleCoreCableWithECCStrategy`?", não uma string específica; qualquer
`mtl_type` futuro que reaproveite essa mesma estratégia herda o dispatch
correto automaticamente. Corrigido também na mesma função que a
`ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios` chama
internamente (mesmo `InternalPerUnitParameters.matrices`), sem precisar de
correção separada em `utils/comsol_data.py` (que já lê `block_sizes` do
resultado, sem checar `mtl_type` — Item 5 de `andreata_case3/BUGS_AND_FIXES.md`).

**Validado:** `andreata_case4.py` roda até o fim sem erro; validação cruzada
contra `andreata_case3` (cenário `'1'`, duto ignorado) com erro relativo
máx. 0,05% (Zs) / 0,12% (Ysh) — mesma ordem de grandeza da validação
equivalente do `andreata_case2` contra `andreata_case1` (0,06%/0,15%).
Regressão limpa em `andreata_case1`, `andreata_case2`, `andreata_case3` e
nos 8 casos irmãos de `underground_flat_model`/`flat_hdpe_enclosed_model`.

**Lição:** ao registrar um novo `mtl_type` que reaproveita uma estratégia
existente (`mtl_strategy_factory`), sempre grepar por checagens de string
exata daquele `mtl_type` em outros módulos (`analytical_forms/`,
`utils/comsol_data.py`, `mtl_main/graphics.py`) — o registro na factory por
si só não garante que todo o pipeline reconheça o novo nome.

---

## Item 3 — Termo próprio do ECC no COMSOL descartado silenciosamente pelo parser de impedância interna

**Onde:** `utils/comsol_data.py` (`ComsolPostProcessor.get_scc_internal_impedance_matrix_combined`).

**Contexto:** ao adicionar `Results/cmsl_internal_impedance_matrix.txt`
(dado COMSOL real, específico do `andreata_case4`, simulação de 3
condutores — núcleo, blindagem, ECC — em vez do arquivo legado de 2
condutores reaproveitado por `andreata_case2`/`hdpe_300mm2`), o gráfico
`internal_impedance_matrix.png` continuou sem nenhum marcador COMSOL para o
elemento próprio do ECC (`p=6, q=6`), mesmo o dado existindo no arquivo.

**Sintoma:** nenhum traceback — só um aviso único no console, `Warning:
COMSOL data not found for 'measured' with key 'impedance_matrix'.`
(`plotter/scc_plotter.py:189`), referente apenas ao componente do ECC (os
três de fase A — `cc`/`cs`/`ss` — plotavam normalmente).

**Causa raiz (duas camadas):**
1. `get_scc_internal_impedance_matrix_combined()` tinha `N = 2` hardcoded e
   lia só 3 das 9 colunas de dado do arquivo novo (`data1`, `data1_1`,
   `data2_1`) — as 6 restantes, incluindo as 3 que envolvem o ECC (mútua
   núcleo-ECC, mútua blindagem-ECC e o termo próprio do ECC, este último na
   última coluna do arquivo, `data3(mf.VCoil_ecc_i0)` → coluna limpa
   `data3_2`), eram descartadas silenciosamente — a matriz retornada
   permanecia `(freq, 2, 2)` independentemente do arquivo ter 4 ou 9
   colunas de dado.
2. Mesmo lendo essas colunas, uma segunda incompatibilidade: o plotter usa
   um único par `(p, q)` para indexar as três fontes de um mesmo componente
   (analítica/COMSOL/MATLAB — `plotter/scc_plotter.py:170-202`), e o
   `plot_config.py` do case4 pede `p=6, q=6` (índice **global** do ECC no
   espaço de 7 condutores, mesma convenção do MATLAB/analítico). Uma matriz
   COMSOL "local" de 3 condutores (núcleo=0, blindagem=1, ECC=2) teria
   `Zi[:, 6, 6]` fora dos limites de qualquer forma — diferente de
   `cc`/`cs`/`ss`, que só funcionam porque a fase A coincide por acaso com
   os índices globais 0/1.

**Correção:** a função agora detecta o formato do arquivo pela presença da
coluna `data1_2` (só existe quando há 3 colunas `data1(...)`, ou seja,
quando o ECC está presente) e, nesse caso, monta a matriz já no espaço
**global** de 7 condutores (`N = 7`, núcleo=0, blindagem=1, ECC=6 — mesma
convenção de `MatlabDataReader.conductor_order` e do `p=6, q=6` usado em
todo `plot_config.py`), preenchendo só os índices 0/1/6 (o resto permanece
zero) a partir das 9 colunas: `data1_2`/`data2_2` para as mútuas
núcleo-ECC/blindagem-ECC, `data3_2` para o termo próprio do ECC. Para o
arquivo legado de 2 condutores (`data1_2` ausente), o comportamento
permanece **idêntico** ao anterior (`N = 2`, mesmas 3 atribuições).

**Validado:**
- `andreata_case4.py`: o aviso de dado COMSOL não encontrado desaparece; o
  gráfico `internal_impedance_matrix.png` passa a mostrar um 4º marcador
  COMSOL (verde, ECC) em `p=6, q=6`, que acompanha de perto os pontos
  MATLAB do ECC (também verdes) em ambos os painéis (R e L) — consistente
  com a mesma ordem de grandeza e formato de curva, mesmo padrão de
  concordância já observado para `cc`/`cs`/`ss`.
- Regressão limpa (mesmo arquivo/função, formato de 2 condutores) em
  `andreata_case1`, `andreata_case2`, `andreata_case3`, `hdpe_300mm2`,
  `scc_132kV_xue`, `scc_34kV_andreata`, `scc_138kV_prysmian` — todos os
  outros consumidores de `load_scc_earth_return_and_internal_scenarios`/
  `get_scc_internal_impedance_matrix_combined` no repositório.

**Lição:** um parser com dimensão de matriz hardcoded (`N = 2`) não avisa
quando o arquivo de entrada cresce (mais colunas) — ele só lê o que já
esperava e descarta o resto em silêncio. Vale detectar o formato pelo
conteúdo (presença/ausência de uma coluna-chave) em vez de assumir um
tamanho fixo, especialmente em parsers de dados externos (COMSOL/MATLAB)
que têm mais de um layout de arquivo em uso no repositório.
