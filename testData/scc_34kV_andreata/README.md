# scc_34kV_andreata

## Propósito

Script de simulação e análise detalhada para um **único cabo coaxial de núcleo único (SCC)** enterrado, com foco na decomposição dos parâmetros PUL em suas contribuições individuais: interna, de retorno à terra e a composição resultante. Calcula impedância série, admitância shunt e coeficientes de potencial por unidade de comprimento ao longo de uma faixa de frequência ampla (1 Hz a 10 MHz), comparando três formulações analíticas.

Diferencia-se do caso `scc_flat_andreata` por analisar **um único cabo isolado**, cobrir uma faixa de frequência mais ampla e apresentar gráficos de decomposição que evidenciam a contribuição relativa de cada parcela dos parâmetros.

---

## Configuração do Sistema (JSON)

Parâmetros carregados de `scc_34kV_andreata.json`:

| Parâmetro | Valor |
|---|---|
| Arranjo | Simples (`single`), 1 cabo |
| Profundidade de enterramento | 1,2 m |
| Solo — condutividade | 0,01 S/m (ρ = 100 Ω·m) |
| Solo — permissividade relativa | 1,0 |
| Núcleo — raio externo | 10,325 mm — condutividade: 38 MS/m |
| Isolação do núcleo | XLPE, 10,675 mm, εr = 2,2 |
| Bainha — raios int./ext. | 21,0 / 21,8 mm — condutividade: 60 MS/m |
| Isolação da bainha | PVC, 2,2 mm, εr = 2,8 |
| Ordem de Fourier | 10 |

O modelo gerado contém **3 condutores**: 1 retorno de solo (`line_id=0`) + 1 par núcleo/bainha (`line_id=1,2`), posicionado em (0, −1,2) m.

---

## Estrutura da Função `main()`

### 1. Construção do modelo MTL (linhas 24–25)

Um único modelo é criado via `SingleCoreCableModelGenerator` + `MulticonductorTransmissionLine`:

| Variável | Solo |
|---|---|
| `mtl_model` | ρ = 100 Ω·m, εr = 1 (referência) |

### 2. Dicionário `pul_data` (linhas 28–48)

Estrutura central que agrega todos os dados da simulação:

```
pul_data
├── 'comsol'            → None (sem dados COMSOL neste caso)
├── 'frequencies'       → 121 pontos log-espaçados em [10⁰, 10⁷] Hz
├── 'internal_matrices' → matrizes internas (preenchido na etapa 3)
└── 'scenarios'         → 3 cenários analíticos (preenchido na etapa 4)
```

**3 cenários analíticos** — mesmo solo, diferentes formulações para retorno à terra:

| Chave | Solo | Formulação Zg | Formulação Yg |
|---|---|---|---|
| `p100_xue` | 100 Ω·m, εr=1 | Magalhães/Xue (integral) | Magalhães/Xue (integral) |
| `p100_deconti` | 100 Ω·m, εr=1 | De Conti (2023) | De Conti (2023) |
| `p100_vance` | 100 Ω·m, εr=1 | De Conti (2023) | Vance (1978) |

### 3. Parâmetros internos analíticos (linhas 50–53)

```python
pul = InternalPerUnitParameters(mtl_model, pul_data['frequencies'])
internal_matrices = pul.matrices()
```

Calcula as matrizes internas do cabo (impedância, admitância e coeficientes de potencial do núcleo e bainha) para as 121 frequências usando a formulação padrão. Resultado armazenado em `pul_data['internal_matrices']`.

> **Diferença em relação ao `scc_flat_andreata`:** a formulação híbrida (`internal_form='hybrid'`) não é usada aqui; o método `matrices()` é chamado com os parâmetros padrão.

### 4. Parâmetros PUL por cenário analítico (linhas 55–63)

Para cada um dos 3 cenários, via `PerUnitParameters`:

1. `earth_return_parameters(zg_form, yg_form)` — calcula impedância **Zg**, coeficiente de potencial **Pg** e admitância de retorno à terra **Yg**.
2. `quasi_tem_approx_matrices(internal_matrices, earth_return)` — monta as matrizes quasi-TEM completas: impedância série **Zs = Zi + Zg**, coeficiente de potencial **Psh = Pi + Pg** e admitância shunt **Ysh**.

Ambos os resultados são armazenados no dicionário de cada cenário.

### 5. Geração de gráficos (linhas 66–85)

`SingleCoreCableModels` gera e salva automaticamente em `Results/` 22 figuras organizadas em três grupos:

#### Grupo 1 — Coeficientes de Potencial

| Método | Arquivo de saída | Conteúdo |
|---|---|---|
| `potential_coefficients_composition('core_sheath')` | `potential_coefficients_composition_core_sheath.png` | Coef. potencial mútuo núcleo-bainha: interno + retorno terra + composição |
| `potential_coefficients_composition('core')` | `potential_coefficients_composition_core.png` | Coef. potencial próprio do núcleo: interno + retorno terra + composição |
| `potential_coefficients_composition('sheath')` | `potential_coefficients_composition_sheath.png` | Coef. potencial próprio da bainha: interno + retorno terra + composição |
| `potential_coefficients_earth_return()` | `earth_return_potential.png` | Coeficiente de potencial de retorno à terra Pg |
| `potential_coefficients_internal()` | `internal_potential.png` | Coeficiente de potencial interno Pi (P00, P01, P11) |

#### Grupo 2 — Impedância Série

| Método | Arquivo de saída | Conteúdo |
|---|---|---|
| `series_impedance_composition('core_sheath')` | `series_impedance_composition_core_sheath.png` | Impedância mútua núcleo-bainha: Zi + Zg + Zs + composição |
| `series_impedance_composition('core')` | `series_impedance_composition_core.png` | Impedância própria do núcleo: Zi + Zg + Zs + composição |
| `series_impedance_composition('sheath')` | `series_impedance_composition_sheath.png` | Impedância própria da bainha: Zi + Zg + Zs + composição |
| `series_impedance_earth_return()` | `earth_return_impedance.png` | Impedância de retorno à terra Zg (elemento [0,0]) |
| `series_impedance_internal()` | `internal_impedance.png` | Impedância interna Zi (própria e mútua núcleo/bainha) |
| `series_impedance_matrix()` | `series_impedance_matrix.png` | Matriz Zs completa (própria núcleo, própria bainha, mútua) |

#### Grupo 3 — Admitância Shunt

| Método | Arquivo de saída | Conteúdo |
|---|---|---|
| `shunt_admittance_composition('core_sheath')` | `shunt_admittance_composition_core_sheath.png` | Admitância mútua núcleo-bainha: Yi + Yg + Ysh + composição |
| `shunt_admittance_composition('core')` | `shunt_admittance_composition_core.png` | Admitância própria do núcleo: Yi + Yg + Ysh + composição |
| `shunt_admittance_composition('sheath')` | `shunt_admittance_composition_sheath.png` | Admitância própria da bainha: Yi + Yg + Ysh + composição |
| `shunt_admittance_earth_return()` | `earth_return_admittance.png` | Admitância de retorno à terra Yg (elemento [0,0]) |
| `shunt_admittance_internal()` | `internal_admittance.png` | Admitância interna Yi (própria e mútua núcleo/bainha) |
| `shunt_admittance_matrix()` | `shunt_admittance_matrix.png` | Matriz Ysh completa (própria núcleo, própria bainha, mútua) |

#### Esquema do sistema

| Método | Arquivo de saída | Conteúdo |
|---|---|---|
| `GroundReturnMTLRepresentation(...)` | `system_schematic.png` | Esquema circuital do sistema MTL (unidades em cm) |

---

## Lógica de Decomposição

Os gráficos de composição sobrepõem quatro curvas para cada elemento da matriz:

| Curva | Cor | Significado |
|---|---|---|
| **Internal** | Azul, tracejado | Contribuição interna do cabo (Zi ou Yi ou Pi) |
| **Earth-return** | Verde-escuro, tracejado | Contribuição de retorno à terra (Zg ou Yg ou Pg) |
| **Series** (quasi-TEM) | Preto, sólido | Resultado quasi-TEM completo (Zs ou Ysh ou Psh) |
| **Composition** | Vermelho, pontilhado | Soma direta interna + retorno terra (Zi+Zg ou 1/(1/Yi+1/Yg)) |

Esta sobreposição permite verificar visualmente em que faixa de frequência cada contribuição domina.

---

## Diferenças em relação ao `scc_flat_andreata`

| Aspecto | `scc_34kV_andreata` | `scc_flat_andreata` |
|---|---|---|
| Arranjo | 1 cabo simples | 3 cabos em arranjo plano |
| Condutores MTL | 3 (solo + núcleo + bainha) | 7 (solo + 3 núcleos + 3 bainhas) |
| Faixa de frequência | 1 Hz a 10 MHz | 10 kHz a 10 MHz |
| Número de cenários | 3 | 9 |
| Dados COMSOL | Não utilizado | Opcional (quando disponível) |
| Formulação interna | Padrão | Híbrida (`hybrid`) |
| Foco dos gráficos | Decomposição por parcela | Comparação entre cenários de solo |
| Plotter | `SingleCoreCableModels` | `SCCPlotter` |

---

## Dependências

| Módulo | Papel |
|---|---|
| `models.single_core_cable.SingleCoreCableModelGenerator` | Constrói o modelo MTL a partir do JSON |
| `mtl_main.source.MulticonductorTransmissionLine` | Representa o sistema MTL multi-condutor |
| `analytical_forms.single_core_cable.InternalPerUnitParameters` | Matrizes internas do cabo (núcleo + bainha) |
| `analytical_forms.single_core_cable.PerUnitParameters` | Matrizes PUL completas com retorno à terra |
| `plotter.scc_models.SingleCoreCableModels` | Geração e salvamento dos gráficos de decomposição |
| `mtl_main.graphics.GroundReturnMTLRepresentation` | Esquema do circuito MTL |

---

## Execução

```bash
python -m testData.scc_34kV_andreata.scc_34kV_andreata
```

**Saída esperada:** 22 arquivos `.png` salvos em `testData/scc_34kV_andreata/Results/`.

---

## Referências

- **Xue, H.** (2018). *General Formulation and Accurate Evaluation of Earth-Return Parameters for Overhead/Underground Cables*. PhD thesis, École Polytechnique de Montréal.
- **Vance, E. F.** (1978). *Coupling to Shielded Cables*. Wiley-Interscience.
- **De Conti, A. et al.** (2023). Formulação para parâmetros de retorno à terra em cabos subterrâneos.
- **Ametani, A., Ohno, T. & Nagaoka, N.** (2015). *Cable System Transients: Theory, Modeling and Simulation*. Wiley-IEEE Press.
