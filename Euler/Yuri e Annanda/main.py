# -*- coding: utf-8 -*-
from edo_pendulo import *

import click
import numpy as np
import csv

METODOS_CHOICE = ['euler_implicito_linear', 'euler_implicito_nao_linear', 'euler_explicito_nao_linear']


def write_csv(t, p, v):
    with open('resultado.csv', 'wb') as file:
        spamwriter = csv.writer(file, delimiter=',')
        spamwriter.writerow(t)
        spamwriter.writerow(p)
        spamwriter.writerow(v)

@click.command()
@click.option('--metodo', type=click.Choice(METODOS_CHOICE), default=METODOS_CHOICE[0], help=u'Nome do método desejado')
@click.option('--tempo_inicial', type=click.FLOAT, default=0., help='Valor inicial do tempo')
@click.option('--tempo_final', type=click.FLOAT, default=10., help='Valor final do tempo')
@click.option('--n_intervalos', type=click.INT, default=100, help='Quantidade de intervalos')
@click.option('--v_inicial', type=click.FLOAT, default=0., help='Valor da velocidade inicial')
@click.option('--theta_inicial', type=click.FLOAT, default=np.pi/24, help='Valor inicial do angulo theta em radianos')
@click.option('--g', type=click.FLOAT, default=10, help='Valor da gravidade')
@click.option('--l', type=click.FLOAT, default=1, help='Valor do comprimento da corda do pendulo')
def main(metodo, tempo_inicial, tempo_final, n_intervalos, v_inicial, theta_inicial, g, l):
    edo = EDOPendulo(
        tempo_inicial,
        tempo_final,
        n_intervalos,
        v_inicial,
        theta_inicial,
        g,
        l
    )
    if metodo == METODOS_CHOICE[0]:
        t, p, v = edo.euler_implicito_linear()
        write_csv(t, p, v)


if __name__ == '__main__':
    main()