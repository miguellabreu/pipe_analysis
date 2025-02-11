# Análise Estrutural de Tubo

Este repositório contém uma análise estrutural de um tubo usando ANSYS APDL, comparando diferentes configurações:
- Número de elementos (2, 10, 20, 50, 100)
- Não-linearidade geométrica (ON/OFF)
- Carregamentos (Pressão e Força)

## Parâmetros
- **Densidade:** 7000 kg/m³
- **Comprimento:** 30 m
- **Diâmetro externo:** 0.20 m
- **Diâmetro interno:** 0.12 m
- **Módulo de Young:** 100 GPa

## Estrutura
- `src/`: Códigos fonte
- `results/`: Resultados das simulações
- `docs/`: Documentação

## Requisitos
- Python 3.8+
- ANSYS APDL

## Configuração do Ambiente

1. Clone o repositório:
    ```bash
    git clone https://github.com/miguellabreu/pipe-analysis.git
    cd pipe-analysis
    ```

2. Crie e ative o ambiente virtual:
    ```bash
    python -m venv venv
    .\venv\Scripts\activate
    ```

3. Instale as dependências:
    ```bash
    pip install -r requirements.txt
    ```
