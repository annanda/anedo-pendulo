# Trabalho 1 da disciplina Análise Numérica de Equações Diferenciais

## Analisando o movimento de um pêndulo 

### Como rodar o programa: 
```
python main.py [OPTIONS]
```

Exemplo:
```
python main.py --metodo euler_implicito_linear --tempo_inicial 0. --tempo_final 10. --n_intervalos 100 --v_inicial 0. --theta_inicial 80. --g 10. --l 1. 
```

#### Valores padrão 

Você pode rodar o programa sem passar nenhuma opção como parâmetro, nesse caso, o programa vai processar com os valores padrões que são os listados abaixo:


--metodo euler_implicito_linear

--tempo_inicial 0.

--tempo final 10.

--n_intervalos 100

--v_inicial 0.

--theta_inicial np.pi/24

--g 10

--l 1

### Saída
O processamento vai gerar uma saída em um arquivo CSV cujo nome é 'resultado.csv' que fica no mesmo nível de diretório do arquivo main.py. 

O arquivo CSV tem a seguinte estrutura: 

**Linha 1**: elementos do vetor t (que representa o tempo) separados por vírgula

**Linha 2**: elementos do vetor p (que representa o angulo teta) separados por vígula 

**Linha 3**: elementos do vetor v (que representa a velocidade) separados por vírgula

### Help
Para obter ajuda sobre os parâmetros e seus tipos, use: 

```
 python main.py --help
```

``` 
Usage: main.py [OPTIONS]

Options:
  --metodo [euler_implicito_linear|euler_implicito_nao_linear|euler_explicito_nao_linear]
                                  Nome do método desejado
  --tempo_inicial FLOAT           Valor inicial do tempo
  --tempo_final FLOAT             Valor final do tempo
  --n_intervalos INTEGER          Quantidade de intervalos
  --v_inicial FLOAT               Valor da velocidade inicial
  --theta_inicial FLOAT           Valor inicial do angulo theta em radianos
  --g FLOAT                       Valor da gravidade
  --l FLOAT                       Valor do comprimento da corda do pendulo
  --help                          Show this message and exit.

```