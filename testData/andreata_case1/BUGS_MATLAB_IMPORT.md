# Resumo de bugs — importação de dados `.mat` (validação `andreata_case1`)

**Contexto:** `testData/andreata_case1/andreata_case1.py` compara resultados analíticos do pyLCP (`InternalPerUnitParameters`) com uma referência exportada do MATLAB (`Results/andreata_internal_impedance_matrix.mat`, variável `Z`, shape `(6,6,90)` — impedância interna complexa de um sistema de 3 cabos unipolares, núcleo+blindagem, arranjo plano). Foram encontrados **4 bugs distintos**, todos na camada de importação/alinhamento dos dados do MATLAB — os dados brutos do `.mat` em si estão corretos.

---

## Bug 1 — Chaves de arquivo `.mat` inexistentes (`matlab_reader.data.get`)

**Onde:** `andreata_case1.py` (antes da correção)

**O quê:** O código buscava `matlab_reader.data.get('andreata_series_impedance_matrix')` e `matlab_reader.data.get('andreata_shunt_admittance_matrix')`. `MatlabDataReader` (`utils/matlab_data.py:63-75`) indexa os dados pelo **nome do arquivo sem extensão** (glob de `Results/*.mat`). Só existe `andreata_internal_impedance_matrix.mat` no diretório — essas duas chaves nunca existiram (resquício de um esquema de nomenclatura anterior). Como o acesso é via `.get()`, o erro é silencioso: retorna `None` sem exceção.

**Impacto:** `pul_data['matlab']['scenarios']['measured']` ficava com ambos os campos `None`. Nos gráficos `z11/z12/z22_internal_vs_matlab`, o `_plot_scc_internal_vs_matlab` (`plotter/scc_plotter.py:154-159`) tentava indexar `None[:, p, q]`, gerava `TypeError`, capturado silenciosamente — resultado: pontos de referência do MATLAB simplesmente não apareciam no gráfico, só um aviso no console.

**Correção:** usar a chave correta (`andreata_internal_impedance_matrix`) e nomear o campo como o plotter espera por padrão (`internal_impedance_matrix`).

---

## Bug 2 — Eixo de frequência assumido incorretamente

**Onde:** `andreata_case1.py:40` (antes da correção)

**O quê:** O código assumia `np.logspace(0, 7, num=90)` (1 Hz–10 MHz) tanto para o cálculo analítico Python quanto para rotular as 90 amostras do `Z` do MATLAB. A varredura real usada no MATLAB — confirmada pelo arquivo `Results/andreata_frequency_range.mat` (variável `freq1`) fornecido posteriormente — é `np.logspace(-2, 7, num=90)` (**0,01 Hz**–10 MHz): mesmo número de pontos e mesmo teto, mas início 2 décadas mais baixo. Verificação numérica: `freq1` bate com `np.logspace(-2,7,90)` com diferença relativa máxima de 4×10⁻¹⁵.

**Impacto:** como `L = Im(Z)/(2π·f)`, dividir cada amostra pelo `f` errado distorce sistematicamente a curva de indutância. Com o eixo errado, `L(f)` aparecia **crescendo monotonicamente ~73×** ao longo da varredura (2,2 nH → 162 nH) — fisicamente impossível para uma auto-impedância interna (deveria ser um platô em baixa frequência seguido de queda por efeito pelicular). Com o eixo correto, a curva mostrou exatamente o comportamento esperado: platô constante (~220 nH) enquanto R ainda está em regime DC, depois decaimento suave até ~162 nH em alta frequência. A parte real (resistência) não é sensível a esse erro (não depende de divisão por `ω`), por isso o problema só ficava visível na indutância — e não aparecia ao plotar `Z` bruto diretamente no MATLAB (que usa seu próprio eixo nativo correto).

**Observação adicional:** `num=90` era o único caso no repositório inteiro com essa contagem de pontos — todos os casos irmãos (`scc_34kV_andreata`, `scc_132kV_xue`, `scc_138kV_prysmian` etc.) usam `num=121`. Mesma classe de bug de copy-paste já identificada em `ohtl_xue_sec43.py`.

**Correção:** carregar `andreata_frequency_range.mat` e usar esse vetor como `pul_data['matlab']['frequencies']`, em vez de reaproveitar o vetor do lado analítico Python. (Nota: o vetor do lado analítico Python também foi atualizado para `np.logspace(-2, 7, num=90)`, então hoje os dois coincidem — mas o código mantém o carregamento explícito do `.mat` de frequência como fonte de verdade, com fallback.)

---

## Bug 3 — `matlab_matrix_key` inconsistente em `plot_config.py`

**Onde:** `testData/andreata_case1/plot_config.py:156` e `:172` (antes da correção)

**O quê:** Os três gráficos `z11/z12/z22_internal_vs_matlab` têm o suptitle "internal-only analytical (Zi) vs. MATLAB reference (Z)", mas apenas `z11` usava `'matlab_matrix_key': 'internal_impedance_matrix'`. `z12` e `z22` ainda apontavam para `'series_impedance_matrix'` — chave que não existe em `pul_data['matlab']['scenarios']['measured']` (consequência direta do Bug 1).

**Impacto:** mesmo depois de corrigir o Bug 1, os gráficos `z12` e `z22` continuariam sem mostrar a referência do MATLAB (mesmo sintoma do Bug 1: `KeyError`/`TypeError` capturado silenciosamente).

**Correção:** alinhadas as três entradas para `'matlab_matrix_key': 'internal_impedance_matrix'`.

---

## Bug 4 — Convenção de ordenação de condutores incompatível entre pyLCP e MATLAB

**Onde:** estrutural — descoberto ao revisar `z12_internal_vs_matlab` e `z22_internal_vs_matlab` após corrigir os Bugs 1–3.

**O quê:** o MATLAB exporta os 6 condutores **agrupados por tipo**: `[core_A, core_B, core_C, sheath_A, sheath_B, sheath_C]` (índices 0–2 = núcleos, 3–5 = blindagens). O pyLCP monta sua matriz interna **agrupada por cabo**: `[core_A, sheath_A, core_B, sheath_B, core_C, sheath_C]`, via `np.kron(np.identity(N), Zij)` em `InternalPerUnitParameters.matrices()` (`analytical_forms/single_core_cable.py:848`), refletindo a ordem de atribuição de `conductor_id` em `models/single_core_cable.py:127-141` (núcleo, depois blindagem, por cabo).

**Evidência:** os elementos diagonais do `Z` do MATLAB confirmam isso — `Z[0,0]==Z[1,1]==Z[2,2]` (≈7,8536e-5 Ω em DC) e `Z[3,3]==Z[4,4]==Z[5,5]` (≈1,549406e-4 Ω em DC), batendo exatamente com a resistência DC teórica (`1/(σ·A)`) do núcleo e da blindagem, respectivamente, calculada a partir do JSON do caso.

**Impacto:** como o código usava o **mesmo par `(p,q)`** para indexar tanto a matriz interna do Python quanto a matriz do MATLAB:
- `z11` (p=0,q=0) coincidia por acaso — índice 0 é `core_A` nas duas convenções.
- `z22` (p=1,q=1): pyLCP pegava a auto-impedância da **blindagem A**; o MATLAB (sem correção) entregava a auto-impedância do **núcleo B** (idêntica ao núcleo A) — grandezas físicas diferentes sendo comparadas lado a lado no mesmo gráfico.
- `z12` (p=0,q=1): pyLCP pegava o acoplamento mútuo **núcleo A ↔ blindagem A** (mesmo cabo, acoplamento coaxial forte); o MATLAB entregava **núcleo A ↔ núcleo B** (cabos diferentes, acoplamento externo fraco) — ordem de grandeza de L incompatível (chegava a ser ~11× maior no par errado).

Validado que os índices corretos (`Z[3,3]` para blindagem A, `Z[0,3]` para núcleo A↔blindagem A) produzem curvas fisicamente consistentes, inclusive convergindo `L(núcleo↔blindagem) → L(blindagem, auto)` em alta frequência — o limite esperado de acoplamento unitário quando o efeito pelicular concentra a corrente nas superfícies enfrentadas.

**Correção:** em `andreata_case1.py:127-130`, a matriz do MATLAB é reordenada com a permutação `[0,3,1,4,2,5]` (aplicada às duas dimensões de condutor) logo no carregamento, antes de ser armazenada em `pul_data['matlab']`, alinhando-a à convenção do pyLCP de uma vez por todas — sem precisar de índices `(p,q)` diferentes por gráfico.

---

## Recomendações para o fluxo de exportação MATLAB → pyLCP

1. **Sempre exportar o vetor de frequência junto com os dados** (como já foi feito para este caso com `andreata_frequency_range.mat`) — evita reincidência do Bug 2 em novos casos de validação.
2. **Documentar/padronizar a convenção de ordenação de condutores** usada nas exportações do MATLAB (por tipo vs. por cabo), ou já exportar na mesma ordem que o pyLCP usa internamente (agrupado por cabo), eliminando a necessidade da permutação de reindexação.
3. Os Bugs 1 e 3 (nomes de chave/arquivo) sugerem que pode valer a pena um teste de sanidade automático no `MatlabDataReader` ou em `andreata_case1.py` que avise (ou falhe alto, em vez de retornar `None` silenciosamente) quando uma chave esperada não é encontrada — o uso de `.get()` sem verificação mascarou os 3 primeiros bugs por bastante tempo.
