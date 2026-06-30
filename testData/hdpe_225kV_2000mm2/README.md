# hdpe_225kV_2000mm2.py — Documentação Completa

## Visão Geral

Script de simulação eletromagnética para um cabo de energia de alta tensão (225 kV, 2000 mm²) enterrado e envolto em um tubo de HDPE (High-Density Polyethylene). O objetivo central é calcular os parâmetros por unidade de comprimento (PUL) — impedância interna e admitância de derivação — por meio de formulações analíticas e compará-los com resultados de referência obtidos no COMSOL Multiphysics.

---

## Problema Físico

Um cabo monopolar (single-core cable, SCC) é enterrado no solo envolto em um tubo de HDPE. A geometria é excêntrica: o eixo do cabo não coincide com o eixo do tubo. O sistema possui:

- **Condutor de fase (core):** condutor de cobre central (área 2000 mm²).
- **Isolação principal:** dielétrico entre core e blindagem metálica.
- **Blindagem metálica (sheath):** condutor tubular externo.
- **Tubo de HDPE:** envoltório de proteção mecânica com folga de ar entre a blindagem e a parede interna do tubo.
- **Solo:** meio condutor de retorno (modelado como plano de terra).

A presença do HDPE cria uma região mista (folga de ar + HDPE) que precisa ser reduzida a um dielétrico equivalente para possibilitar o uso das formulações analíticas padrão do tipo MTL (Multiconductor Transmission Line).

---

## Cenários de Modelagem

O script define quatro modelos MTL (`model_0` a `model_3`) que representam diferentes hipóteses sobre como tratar a região entre a blindagem e o solo:

| Modelo | Variável | Descrição |
|--------|----------|-----------|
| `model_0` | `mtl_0` | Cabo excêntrico com tubo HDPE — geometria real completa. Usado como referência geométrica para o cálculo dos ERS. |
| `model_1` | `mtl_1` | Cabo enterrado sem tubo HDPE — ignora a presença do envoltório. Usado como Cenário 1 (underground). |
| `model_2` | `mtl_2` | Cabo sem enclosure; a região entre blindagem e solo é substituída por um dielétrico equivalente calculado pelo método de **ponderação por área** (area-weighted `epsr`). |
| `model_3` | `mtl_3` | Cabo sem enclosure; permissividade equivalente calculada pelo método **GMD** (Geometric Mean Distance), caso 3.1. |

Os Cenários 1, 2 e 3 correspondem respectivamente a `mtl_1`, `mtl_2` e `mtl_3` na estrutura `pul_data['scenarios']`.

---

## Faixa de Frequência

```python
frequencies = np.logspace(0, 6, num=31)  # 1 Hz a 1 MHz, 31 pontos em escala logarítmica
```

---

## Dependências

### Módulos da Biblioteca Padrão

| Módulo | Uso |
|--------|-----|
| `sys` | Encerramento em caso de erro de importação |
| `os` | Limpeza do terminal (`cls`/`clear`) |
| `copy` | Cópia profunda de dicionários de modelo |
| `time` | Medição do tempo de execução |
| `numpy` | Operações numéricas e geração de frequências |
| `matplotlib.pyplot` | Exibição dos gráficos ao final |

### Módulos PyLCP

| Módulo | Classe / Símbolo | Responsabilidade |
|--------|------------------|-----------------|
| `utils.case_utils` | `*` (wildcard) | Utilitários de formatação, impressão de matrizes e salvamento de figuras |
| `utils.comsol_data` | `ComsolPostProcessor` | Leitura e processamento dos resultados exportados do COMSOL |
| `plotter.scc_plotter` | `HDPEPlotter` | Geração dos gráficos de impedância interna |
| `mtl_main.graphics` | `GroundReturnMTLRepresentation` | Geração dos esquemáticos de seção transversal |
| `mtl_main.source` | `MulticonductorTransmissionLine` | Construção e validação do modelo MTL |
| `models.single_core_cable` | `SingleCoreCableModelGenerator` | Geração dos dicionários de modelo a partir do JSON do caso |
| `analytical_forms.single_core_cable` | `InternalPerUnitParameters` | Cálculo analítico das matrizes PUL (impedância e admitância) |
| `analytical_forms.single_core_cable` | `EquivalentRadiiSystems` | Cálculo dos raios e permissividades equivalentes |
| `.plot_config` | `PLOT_CONFIG` | Configuração dos gráficos (séries, estilos, eixos) |

---

## Estrutura de Dados Principal — `pul_data`

```python
pul_data = {
    'frequencies': np.ndarray,          # (31,) — vetor de frequências

    'comsol': {
        'frequencies': np.ndarray,       # frequências do arquivo COMSOL
        'angular_frequencies': np.ndarray,
        'scenarios': {
            '1': {
                'coaxial_cable_impedance': dict,       # z11, z12, z2i, Zcs
                'internal_impedance_matrix': dict,     # js_method, energy_method
                'internal_impedance_elements': dict,   # self_core, self_sheath, mutual
            }
        }
    },

    'scenarios': {
        '1': {
            'mtl': MulticonductorTransmissionLine,     # modelo underground
            'internal_parameters': dict,               # saída de parameters_hybrid()
            'internal_matrices': dict,                 # saída de matrices()
        },
        '2': { ... },   # modelo area-weighted ERS
        '3': { ... },   # modelo GMD ERS
    },

    'analytical': {
        'frequencies': np.ndarray,
        'scenarios': { ... }            # referência ao mesmo dict de 'scenarios'
    }
}
```

---

## Fluxo de Execução — `main()`

### 1. Geração dos Modelos Geométricos

```python
model_generator = SingleCoreCableModelGenerator(__file__)
model_0 = model_generator.eccentric_hdpe_enclosed_model()
model_1 = model_generator.underground_model()
```

- `SingleCoreCableModelGenerator` carrega os parâmetros físicos do arquivo `hdpe_225kV_2000mm2.json` no mesmo diretório.
- `eccentric_hdpe_enclosed_model()` retorna um dicionário com condutores, isolações e enclosure (tubo HDPE) posicionados excentricamente.
- `underground_model()` retorna o mesmo cabo sem enclosure, representando o enterramento simples.

### 2. Construção dos Objetos MTL

```python
mtl_0 = MulticonductorTransmissionLine(model_0)
mtl_1 = MulticonductorTransmissionLine(model_1)
```

`MulticonductorTransmissionLine` valida o modelo, extrai superfícies condutoras e calcula propriedades geométricas usadas pelas formulações analíticas.

### 3. Cálculo dos Sistemas de Raios Equivalentes (ERS)

#### Método Area-Weighted → `model_2`

```python
ers = EquivalentRadiiSystems(mtl_0)
epsr_area = ers.equiv_rel_permittivity_epsr_area_weighted()
r4 = epsr_area['sheath_outer_radius']
r7 = epsr_area['sheath_enclosure_outer_radius']

model_2 = copy.deepcopy(model_0)
model_2[2]['enclosure'] = None
model_2[2]['insulation']['thickness'] = r7 - r4
model_2[2]['insulation']['relative_permittivity'] = epsr_area['equivalent_relative_permittivity']
mtl_2 = MulticonductorTransmissionLine(model_2)
```

A permissividade equivalente é calculada por média ponderada pela área das regiões dielétrico, folga de ar e HDPE. O enclosure é removido e substituído por uma camada única de isolação equivalente de espessura `r7 - r4`.

#### Método GMD (Caso 3.1) → `model_3`

```python
ers = EquivalentRadiiSystems(mtl_0)
gmp = ers.equivalent_parameters_from_gmd()
eps_a = gmp['equivalent_relative_permittivity']['case 3.1']

model_3 = copy.deepcopy(model_0)
model_3[2]['enclosure'] = None
model_3[2]['insulation']['relative_permittivity'] = eps_a
mtl_3 = MulticonductorTransmissionLine(model_3)
```

Baseado na metodologia GMD de Lafaia (2015). No caso 3.1, o raio externo do isolador equivalente é o raio externo do cabo (`r5`), sem modificar a espessura da camada.

### 4. Leitura dos Dados COMSOL

```python
cmsl_processor = ComsolPostProcessor(__file__)
cmsl_params = cmsl_processor.get_general_parameters('cmsl_coaxial_cable_impedance')
pul_data['comsol'].update(cmsl_params)

value['coaxial_cable_impedance'] = cmsl_processor.get_coaxial_cable_parameters()
scc_elements = cmsl_processor.get_scc_internal_impedance_elements()
value['internal_impedance_matrix'] = scc_elements
value['internal_impedance_elements'] = scc_elements
```

`ComsolPostProcessor` localiza automaticamente o diretório `Results/` do caso e lê os arquivos exportados. Os dados obtidos são:

| Chave | Conteúdo |
|-------|----------|
| `coaxial_cable_impedance` | `z11`, `z12`, `z2i`, `Zcs` — elementos de impedância de cabo coaxial |
| `internal_impedance_matrix` | Matrizes 2×2 de impedância interna (método JS e método da energia) |
| `internal_impedance_elements` | Elementos individuais: `self_core`, `self_sheath`, `mutual` |

### 5. Cálculo Analítico dos Parâmetros PUL

```python
for key, value in pul_data['scenarios'].items():
    pul = InternalPerUnitParameters(value['mtl'], pul_data['frequencies'])
    value['internal_parameters'] = pul.parameters_hybrid()
    value['internal_matrices'] = pul.matrices()
```

Para cada cenário analítico (1, 2, 3):

- `parameters_hybrid()` — abordagem híbrida: funções de Bessel abaixo de 100 kHz, aproximações acima. Retorna os elementos individuais `z11`, `z12`, `z2i`, `Zcs`, etc.
- `matrices()` — monta as matrizes completas de impedância série e admitância shunt no domínio da frequência.

**Estrutura de retorno de `parameters_hybrid()`:**

```python
{
    'zcs': {'z11': array, 'z12': array, 'z2i': array, 'Zcs': array},
    'zs3': {'z20': array, 'z23': array},
    'z2m': array,
    'potentials': {'pcj': array, 'psj': array}
}
```

**Estrutura de retorno de `matrices()`:**

```python
{
    'impedance_matrix':          np.ndarray,  # (31, N, N) — Ω/m
    'resistance_matrix':         np.ndarray,  # (31, N, N) — Re(Z)
    'inductance_matrix':         np.ndarray,  # (31, N, N) — Im(Z)/(2πf)
    'shunt_admittance_matrix':   np.ndarray,  # (31, N, N) — S/m
    'potential_coefficient_matrix': np.ndarray,  # (N, N) — independente da frequência
    'capacitance_matrix':        np.ndarray,  # (N, N) — F/m
}
```

### 6. Geração dos Gráficos

#### Gráficos de Impedância Interna

```python
plotter = HDPEPlotter(__file__, pul_data, PLOT_CONFIG, autoSave=False)
plotter.hdpe_internal_impedance_matrix()
plotter.hdpe_internal_impedance_elements()
```

- `hdpe_internal_impedance_matrix()` — gráfico duplo (resistência e indutância) comparando os três cenários analíticos com os dados COMSOL para os elementos da matriz PUL (Z_cc, Z_cs, Z_ss).
- `hdpe_internal_impedance_elements()` — decompõe a impedância em contribuições individuais de cada camada (z11, z12, z2i), comparando método JS com método da energia.

#### Esquemáticos de Seção Transversal

```python
schematic_configs = [
    {'mtl': mtl_0, 'filename': 'schematic_original_hdpe'},
    {'mtl': mtl_1, 'filename': 'schematic_ignored_hdpe'},
]
for config in schematic_configs:
    schematic = GroundReturnMTLRepresentation(
        __file__, config['mtl'], autoSave=True, units='millimeter'
    )
    schematic.system_schematic(base_filename=config['filename'])
```

Gera dois esquemáticos de seção transversal em `Results/`:
- `schematic_original_hdpe.*` — seção real com o tubo HDPE.
- `schematic_ignored_hdpe.*` — seção sem o tubo (underground).

`autoSave=True` salva os arquivos automaticamente; `autoSave=False` nos plotters de impedância deixa o controle para o `plt.show()` ao final.

---

## Saídas Geradas

| Arquivo | Diretório | Descrição |
|---------|-----------|-----------|
| `schematic_original_hdpe.png/svg` | `Results/` | Seção transversal com tubo HDPE |
| `schematic_ignored_hdpe.png/svg` | `Results/` | Seção transversal sem tubo |
| Gráficos de impedância | Janela interativa | Exibidos via `plt.show()` (não salvos automaticamente) |

---

## Parâmetros Geométricos de Referência

Os raios abaixo identificam as superfícies do sistema em ordem crescente a partir do eixo do cabo:

| Símbolo | Superfície |
|---------|-----------|
| `r1` | Raio interno do core |
| `r2` | Raio externo do core |
| `r3` | Raio externo da isolação do core |
| `r4` | Raio externo da blindagem metálica (sheath) |
| `r5` | Raio externo da isolação da blindagem |
| `r6` | Raio interno do tubo HDPE (início da folga de ar) |
| `r7` | Raio externo do tubo HDPE |

---

## Execução

```bash
# A partir do diretório raiz do projeto (c:\git\PyLCP)
python -m testData.hdpe_225kV_2000mm2.hdpe_225kV_2000mm2
```

O terminal é limpo automaticamente no início de cada execução (`main()`). O tempo total é reportado ao final da simulação.
