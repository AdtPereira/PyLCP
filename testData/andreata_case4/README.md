# andreata_case4

## Propósito

Script de simulação e validação para a **Configuração 4** (Figura 5.4) da
referência: três cabos coaxiais de núcleo único (SCC — núcleo + blindagem)
em arranjo plano, cada um instalado dentro do seu próprio duto de HDPE
(idêntico ao `andreata_case2`), **mais um cabo de aterramento isolado (ECC)
compartilhando o duto do terceiro cabo** — a mesma posição relativa ao cabo
C que o ECC já tem no `andreata_case3` (Configuração 3), só que agora dentro
de um duto. Calcula os parâmetros por unidade de comprimento (PUL) —
internos (núcleo/blindagem/ECC) e de retorno à terra — comparando três
formas analíticas distintas de representar o efeito dielétrico do duto, e
faz validação cruzada com os dados MATLAB do `andreata_case3`.

Referência: Andreata, Luis Eduardo Batista. *Análise das Características de
Propagação e de Transitórios Eletromagnéticos em Cabos Subterrâneos
Instalados em Tubos Não Metálicos no Contexto de Parques Eólicos.* Programa
de Pós-Graduação em Engenharia Elétrica, UFMG, 2025.
https://hdl.handle.net/1843/2086

---

## Configuração do Sistema (JSON)

Parâmetros carregados de `andreata_case4.json`:

| Parâmetro | Valor |
|---|---|
| Arranjo | Plano (`flat`), 3 cabos SCC + 1 ECC |
| Profundidade de enterramento | 1,2 m (centro de cada cabo SCC) |
| Espaçamento entre cabos SCC | 0,2 m (centro a centro) |
| Solo — condutividade | 0,01 S/m (ρ = 100 Ω·m) |
| Solo — permissividade relativa | 1,0 |
| Núcleo — raio externo | 10,325 mm — condutividade: 38 MS/m |
| Isolação do núcleo | XLPE, 10,675 mm, εr = 2,2 |
| Bainha — raios int./ext. | 21,0 / 21,8 mm — condutividade: 60 MS/m |
| Isolação da bainha | PVC, 2,2 mm, εr = 2,8 |
| Duto (HDPE) — raios int./ext. | 50,8 / 57,8 mm — εr = 2,35, excêntrico (cabo apoiado no fundo), **idêntico ao `andreata_case2`** |
| ECC — raio externo | 3,3 mm — condutividade: 38 MS/m, isolação XLPE 2,0 mm |
| Posição do ECC | `ecc_horizontal_gap = 0,02802` m / `ecc_vertical_gap = 0,00886` m (centro-a-centro em relação ao cabo C) — **os mesmos valores de `andreata_case3.json`**, sem recálculo: a posição do ECC em relação ao cabo C não muda entre a Configuração 3 e a 4, só a presença do duto |
| Ordem de Fourier | 10 |

O modelo heterogêneo (usado nos 3 cenários analíticos) contém **7
condutores**: 1 retorno de solo (`line_id=0`) + 3 pares núcleo/bainha
(`line_id=1..6`) + 1 ECC (`line_id=7`) — mesma convenção de IDs do
`andreata_case3`.

---

## Decisão de arquitetura: dois modelos físicos separados

`EquivalentRadiiSystems` (usado para calcular os raios/permissividades
equivalentes ERS/GMD) exige `self.model.scc` como dict **achatado** (uma
seção transversal só) — incompatível com o dict **agrupado por cabo** que a
estratégia heterogênea (`SingleCoreCableWithECCStrategy`, necessária por
causa do ECC) produz. Por isso este caso usa **dois modelos físicos
distintos**, cada um só para o que precisa:

- **`flat_hdpe_enclosed_model()`** (o mesmo método do `andreata_case2`,
  homogêneo, sem ECC) — usado **só** para alimentar `EquivalentRadiiSystems`
  e obter `r4..r7`/permissividades equivalentes. O efeito dielétrico do duto
  não depende da presença do ECC dentro dele, então esse modelo homogêneo é
  suficiente para essa finalidade.
- **`flat_hdpe_enclosed_with_shared_ecc_model()`** (novo método,
  `models/single_core_cable.py`) — modelo físico completo (duto real nas 3
  fases + ECC na posição centro-a-centro herdada do `andreata_case3`) usado
  **só para o esquemático** (`system_schematic.png`). **Nunca** passa por
  `InternalPerUnitParameters`/`PerUnitParameters` — ver limitação abaixo.

Os 3 cenários analíticos (Underground/ERS/GMD) partem, em vez disso, de
`flat_scc_with_ecc_cable_model()` (o mesmo método do `andreata_case3`,
heterogêneo, sem duto), com o efeito do duto aproximado via substituição de
isolação equivalente nas 3 bainhas SCC (`_override_sheath_insulation`,
mesma técnica do `andreata_case2`).

---

## Limitações e Dados Pendentes

- **Não há curva analítica pyLCP para o acoplamento SCC↔ECC dentro do duto
  compartilhado.** A geometria de dois condutores excêntricos (SCC e ECC)
  dividindo um duto HDPE é **não-canônica** — não existe formulação
  analítica fechada para essa configuração específica no método
  quasi-TEM/GMD do pyLCP. Não é uma lacuna de implementação a ser corrigida:
  é um limite físico do método. Por isso os casos-irmãos de cabo único
  `hdpe_ecc_300mm2`/`hdpe_ecc_2000mm2` (que modelam exatamente essa
  geometria compartilhada, via `hdpe_shared_enclosed_model()`) também nunca
  chamam `InternalPerUnitParameters`/`PerUnitParameters` sobre esse modelo —
  só o usam para o esquemático e para comparação direta com dados de
  elementos finitos (COMSOL). Este caso segue o mesmo padrão: a única
  referência real para o efeito do duto sobre o ECC é a comparação direta
  com MATLAB (ver gráficos `*_ecc`).
- **Os 3 cenários (Underground/ERS/GMD) nunca alteram o ECC em si** — a
  substituição de isolação equivalente (`_override_sheath_insulation`) só
  atua nas 3 bainhas SCC. Por isso, nos gráficos `self_impedance_ecc`,
  `self_admittance_ecc`, `earth_return_impedance_ecc` e
  `earth_return_potential_coeff_ecc`, as três curvas analíticas tendem a
  coincidir quase totalmente — comportamento esperado, não um bug.
- **Retorno à terra nunca modela o duto**, em nenhum dos 3 cenários — mesma
  limitação já documentada no `andreata_case2`: a formulação de `Zg`/`Yg`
  enxerga apenas posição/raio externo de cada cabo, nunca a presença do
  duto. A diferença entre os 3 cenários no retorno à terra das fases SCC
  está só no raio externo efetivo usado no termo próprio de imagem
  (estendido até `r7` por ERS/GMD).
- **Dado COMSOL de impedância interna já disponível**
  (`Results/cmsl_internal_impedance_matrix.txt`) — simulação própria do
  case4, 3 condutores (núcleo/blindagem/ECC), incluindo o termo próprio do
  ECC (`p=6, q=6` no gráfico `internal_impedance_matrix.png`). Lido por
  `ComsolPostProcessor.get_scc_internal_impedance_matrix_combined()`, que
  detecta esse formato de 3 condutores (vs. o legado de 2 condutores de
  `andreata_case2`/`hdpe_300mm2`) pela presença da coluna `data1_2` e monta
  a matriz já no espaço global de 7 condutores (núcleo=0, blindagem=1,
  ECC=6) — ver Item 3 de `BUGS_AND_FIXES.md`.
- **Ainda não há dado COMSOL de retorno à terra nem de admitância interna**
  (`cmsl_ground_return_impedance.txt`, `cmsl_internal_admittance_matrix.txt`)
  — o bloco `ComsolPostProcessor.load_scc_earth_return_and_internal_scenarios(...)`
  é chamado mesmo assim (padrão "guardado" do `andreata_case3`), só emite
  avisos de dado ausente para essas duas partes, sem bloquear a execução.
- **MATLAB do duto já disponível.** `Results/andreata_case4_*.mat` — mesmo
  formato `(7,7,90)` tipo-agrupado dos outros casos, `conductor_order=[0, 3,
  1, 4, 2, 5, 6]` (ECC sem par para trocar, fica em último em ambas as
  convenções — mesma ordem do `andreata_case3`).

---

## Execução

```bash
python -m testData.andreata_case4.andreata_case4
```

**Saída esperada:** validação cruzada contra `andreata_case3` impressa no
console (erro relativo máx. ~0,05%/0,12% em `Zs`/`Ysh`, cenário Underground),
`internal_impedance_matrix.png`, `internal_admittance_matrix.png` e
`system_schematic.png` salvos em `testData/andreata_case4/Results/`. Escopo
reduzido deliberadamente por ora (`main()` mantém o restante dos gráficos —
comparação por modelo de duto e retorno à terra — comentado em
`plotter.compare_complete_matrices(...)`, pronto para reativar); o dado
COMSOL de impedância interna já entra nos dois gráficos ativos, o resto será
ligado conforme mais dados COMSOL forem chegando (retorno à terra,
admitância interna).

---

## Referências

- **Andreata, L. E. B.** (2025). *Análise das Características de Propagação e de
  Transitórios Eletromagnéticos em Cabos Subterrâneos Instalados em Tubos Não Metálicos
  no Contexto de Parques Eólicos.* UFMG. https://hdl.handle.net/1843/2086
- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return
  Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de
  Montréal. (formulação `magalhaes_xue` de retorno à terra)
- **Lafaia, I.** (2015). Método GMD para permissividade equivalente de dutos (caso 3.1).
