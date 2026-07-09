# scc_flat_andreata

## Propósito

Script de simulação e validação para um sistema de **três cabos coaxiais de núcleo único (SCC)** enterrados em arranjo plano (*flat*). Calcula os parâmetros por unidade de comprimento (PUL) — impedância série e admitância shunt — para múltiplos cenários de solo e formulações analíticas, e gera os gráficos de referência comparando com resultados COMSOL (quando disponíveis).

O caso reproduz as **Figuras 4.19, 4.21a/b e 4.23, 4.25a/b** da tese de Xue (2018).

---

## Configuração do Sistema (JSON)

Parâmetros carregados de `scc_flat_andreata.json`:

| Parâmetro | Valor |
|---|---|
| Arranjo | Plano (`flat`), 3 cabos |
| Profundidade de enterramento | 1,2 m |
| Espaçamento entre cabos | 0,2 m |
| Solo — condutividade | 0,01 S/m (ρ = 100 Ω·m) |
| Solo — permissividade relativa | 1,0 |
| Núcleo — raio externo | 10,325 mm — condutividade: 38 MS/m |
| Isolação do núcleo | XLPE, 10,675 mm, εr = 2,2 |
| Bainha — raios int./ext. | 21,0 / 21,8 mm — condutividade: 60 MS/m |
| Isolação da bainha | PVC, 2,2 mm, εr = 2,8 |
| Ordem de Fourier | 10 |

O modelo gerado contém **7 condutores**: 1 retorno de solo (`line_id=0`) + 3 pares núcleo/bainha (`line_id=1..6`) posicionados em (0,0), (0,2; 0) e (0,4; 0) m, todos à profundidade −1,2 m.

---

## Estrutura da Função `main()`

### 1. Construção dos modelos MTL (linhas 27–36)

Três variantes do modelo são criadas via `SingleCoreCableModelGenerator` + `MulticonductorTransmissionLine`, representando os cenários de solo:

| Variável | Solo |
|---|---|
| `mtl_model_a` | ρ = 100 Ω·m, εr = 1 (referência) |
| `mtl_model_b` | ρ = 100 Ω·m, εr = 20 |
| `mtl_model_c` | ρ = 500 Ω·m, εr = 1 |

### 2. Dicionário `pul_data` (linhas 38–94)

Estrutura central que agrega todos os dados da simulação:

```
pul_data
├── 'frequencies'       → 121 pontos log-espaçados em [10⁴, 10⁷] Hz
├── 'internal_matrices' → matrizes internas (preenchido na etapa 4)
├── 'comsol'            → dados COMSOL (preenchido se arquivo disponível)
│   └── 'scenarios'     → 3 cenários COMSOL (rho_g_100/500, epsr1_1/20)
└── 'scenarios'         → 9 cenários analíticos (preenchido na etapa 5)
```

**9 cenários analíticos** — combinação de 3 solos × 3 formulações:

| Chave | Solo | Formulação Zg | Formulação Yg |
|---|---|---|---|
| `p100_er1` | 100 Ω·m, εr=1 | Magalhães/Xue | Magalhães/Xue |
| `p100_er20` | 100 Ω·m, εr=20 | Magalhães/Xue | Magalhães/Xue |
| `p500_er1` | 500 Ω·m, εr=1 | Magalhães/Xue | Magalhães/Xue |
| `p100_er1_vance` | 100 Ω·m, εr=1 | De Conti | Vance (1978) |
| `p100_er20_vance` | 100 Ω·m, εr=20 | De Conti | Vance (1978) |
| `p500_er1_vance` | 500 Ω·m, εr=1 | De Conti | Vance (1978) |
| `p100_er1_deconti` | 100 Ω·m, εr=1 | De Conti | De Conti (2023) |
| `p100_er20_deconti` | 100 Ω·m, εr=20 | De Conti | De Conti (2023) |
| `p500_er1_deconti` | 500 Ω·m, εr=1 | De Conti | De Conti (2023) |

### 3. Processamento COMSOL — opcional (linhas 96–114)

Tenta carregar `Results/cmsl_ground_return_impedance.txt` via `ComsolPostProcessor`. Se o arquivo existir, constrói as matrizes internas nas frequências COMSOL e calcula os parâmetros de retorno à terra para cada um dos 3 cenários COMSOL. Se não existir, o bloco é ignorado com aviso e `pul_data['comsol']` é esvaziado.

### 4. Parâmetros internos analíticos (linhas 116–119)

```python
pul = InternalPerUnitParameters(mtl_model_a, pul_data['frequencies'])
internal_matrices = pul.matrices(internal_form='hybrid')
```

Calcula, para as 121 frequências, as matrizes internas do cabo (impedância do núcleo, bainha e acoplamentos) usando a formulação híbrida. Resultado armazenado em `pul_data['internal_matrices']`.

### 5. Parâmetros PUL por cenário analítico (linhas 121–128)

Para cada um dos 9 cenários, via `PerUnitParameters`:

1. `earth_return_parameters(zg_form, yg_form)` — calcula a impedância e admitância de retorno à terra **Zg** e **Yg** usando a formulação especificada.
2. `quasi_tem_approx_matrices(internal_matrices, earth_return)` — monta as matrizes quasi-TEM completas: impedância série **Zs = Zi + Zg** e admitância shunt **Ysh**.

Ambos os resultados são armazenados no dicionário do cenário.

### 6. Geração de gráficos (linhas 131–138)

`SCCPlotter` gera e salva automaticamente em `Results/` os seguintes gráficos:

| Método | Figuras | Conteúdo |
|---|---|---|
| `scc_earth_propagation_constant()` | `earth_propagation_constant.png` | Constante de propagação do solo γ₁ (α e β) |
| `scc_earth_return_impedance_matrix()` | `earth_return_impedance_self/mutual_ab/ac.png` | Zg — autoimpedância e mútuos AB, AC |
| `scc_earth_return_admittance_matrix()` | `earth_return_admittance_self/mutual_ab/ac.png` | Yg — autoadmitância e mútuos AB, AC |
| `scc_series_impedance_matrix(...)` | `fig419`, `fig421a`, `fig421b` | Zs — autoimpedância e mútuos da bainha |
| `scc_shunt_admittance_matrix(...)` | `fig423`, `fig425a`, `fig425b` | Ysh — autoadmitância e mútuos da bainha |
| `GroundReturnMTLRepresentation(...)` | `system_schematic.png` | Esquema circuital do sistema MTL |

---

## Dependências

| Módulo | Papel |
|---|---|
| `models.single_core_cable.SingleCoreCableModelGenerator` | Constrói o modelo MTL a partir do JSON |
| `mtl_main.source.MulticonductorTransmissionLine` | Representa o sistema MTL multi-condutor |
| `analytical_forms.single_core_cable.InternalPerUnitParameters` | Matrizes internas do cabo (núcleo + bainha) |
| `analytical_forms.single_core_cable.PerUnitParameters` | Matrizes PUL completas com retorno à terra |
| `utils.comsol_data.ComsolPostProcessor` | Leitura e processamento dos dados COMSOL |
| `plotter.scc_plotter.SCCPlotter` | Geração e salvamento dos gráficos |
| `mtl_main.graphics.GroundReturnMTLRepresentation` | Esquema do circuito MTL |
| `plot_config.PLOT_CONFIG` | Layout, séries e limites de todos os gráficos |

---

## Execução

```bash
python -m testData.scc_flat_andreata.scc_flat_andreata
```

**Saída esperada:** 11 arquivos `.png` salvos em `testData/scc_flat_andreata/Results/`. Os dados COMSOL são opcionais — a simulação analítica roda normalmente sem eles.

---

## Referências

- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de Montréal.
- **Vance, E. F.** (1978). *Coupling to Shielded Cables*. Wiley-Interscience.
- **De Conti, A. et al.** (2023). Formulação para parâmetros de retorno à terra em cabos subterrâneos.
