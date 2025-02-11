
from ansys.mapdl.core import launch_mapdl
import numpy as np
import os
import csv

# Parâmetros
rho_s = 7000.0      # Densidade do riser em kg/m³
L_s = 30.0          # Comprimento do riser em metros
De_s = 0.20         # Diâmetro externo do riser em metros
Di_s = 0.12         # Diâmetro interno do riser em metros
E_s = 100.0e9       # Módulo de Young do riser em Pa

Thic_riser = (De_s - Di_s) / 2  # Espessura do riser

pi = np.pi
A_s = 0.25 * pi * (De_s**2 - Di_s**2)  # Área da seção transversal
grav_z = 9.81                          # Aceleração gravitacional em m/s²
Force1 = rho_s * A_s * L_s * grav_z    # Força total
q = Force1 / L_s                       # Carga distribuída
Inercia = pi * (De_s**4 - Di_s**4) / 64  # Momento de inércia
dzin = (5 * q * (L_s**4)) / (384 * E_s * Inercia)  # Deslocamento teórico máximo

# Diretório base de saída
base_dir = "C:/Users/Miguel/Desktop/Lamce/pipe/force"

# Casos a serem executados
nels = [2, 10, 20, 50, 100]
nlgeom_options = ['OFF', 'ON']

# Loop sobre nel e nlgeom_options
for nel in nels:
    for nlgeom_option in nlgeom_options:
        # Diretório único para esta execução
        run_dir = os.path.join(base_dir, f'nel_{nel}_nlgeom_{nlgeom_option}')
        os.makedirs(run_dir, exist_ok=True)

        # Inicializar o MAPDL dentro do loop
        mapdl = launch_mapdl(run_location=run_dir, override=True)
       
        # Início da modelagem
        mapdl.clear()
        mapdl.prep7()
        # Definição da geometria
        mapdl.k(1, 0, 0, 0)
        mapdl.k(2, L_s, 0, 0)
        mapdl.l(1, 2)
        # Propriedades do material
        matpipe = 1
        mapdl.mp('EX', matpipe, E_s)
        mapdl.mp('PRXY', matpipe, 0.3)
        mapdl.mp('DENS', matpipe, rho_s)
        # Definição do elemento PIPE288
        mapdl.et(1, 'PIPE288')
        mapdl.sectype(1, 'PIPE', 'riser')
        mapdl.secdata(De_s, Thic_riser, 12)
        mapdl.keyopt(1, 1, 1)  # Pipe espesso
        # Malha
        Esize_riser = L_s / nel
        mapdl.lesize(1, Esize_riser)
        mapdl.mshkey(1)
        mapdl.lmesh('ALL')
        # Condições de contorno
        mapdl.d(1, 'UX', 0)  # Nó inicial fixo
        mapdl.d(1, 'UY', 0)
        mapdl.d(1, 'UZ', 0)
        mapdl.d(1, 'ROTX', 0)
        mapdl.d(2, 'UX', 0)  # Nó final fixo
        mapdl.d(2, 'UY', 0)
        mapdl.d(2, 'UZ', 0)
        # Salvar o modelo
        model_name = f'riser_model_nel_{nel}_nlgeom_{nlgeom_option}'
        mapdl.save(model_name)
        mapdl.finish()
        print(f"Modelo salvo como '{model_name}.db'.")
        # Solução
        mapdl.slashsolu()
        mapdl.antype('STATIC')
        mapdl.nlgeom(nlgeom_option)
        
        # Aplicar carga concentrada
        
        fj = -(Force1 / L_s) * Esize_riser  # Carga concentrada em cada nó

        mapdl.f('ALL', 'Fy', fj)  # Aplica a carga concentrada no nó selecionado (na direção Y)

        mapdl.solve()
        mapdl.finish()
        # Pós-processamento
        mapdl.post1()
        mapdl.set(1, 1)
         # Obter deslocamentos nodais
        displacements = mapdl.post_processing.nodal_displacement('ALL')
        # Salvar deslocamentos em CSV
        csv_filename = f'displacements_nel_{nel}_nlgeom_{nlgeom_option}.csv'
        csv_filepath = os.path.join(run_dir, csv_filename)
        with open(csv_filepath, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Node', 'UX', 'UY', 'UZ'])
            for displacement in displacements:
                node = displacement[0]
                if len(displacement) >= 4:
                    ux, uy, uz = displacement[1:4]
                else:
                    ux, uy, uz = displacement[1], 0, 0  # Ajuste conforme necessário
                writer.writerow([int(node), ux, uy, uz])
        print(f"Deslocamentos salvos em '{csv_filepath}'.")
        # Obter deslocamento no nó central
        mapdl.nsel('S', 'LOC', 'X', L_s / 2)
        noc = mapdl.get(entity='NODE', item1='NUM', it1num='MAX')
        disp_x = mapdl.get_value('NODE', noc, 'U', 'X')
        disp_y = mapdl.get_value('NODE', noc, 'U', 'Y')
        disp_z = mapdl.get_value('NODE', noc, 'U', 'Z')
        print(f"Deslocamento no nó central (Nó {noc}): UX = {disp_x}, UY = {disp_y}, UZ = {disp_z}")
        # Adicionar deslocamento do nó central ao CSV
        with open(csv_filepath, mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Nó Central', disp_x, disp_y, disp_z])
        # Salvar o modelo após a solução
        mapdl.save(f'{model_name}_after_solve')
        mapdl.finish()
        print(f"Modelo salvo como '{model_name}_after_solve.db'.")
        # Instruções para abrir no ANSYS MAPDL
        print("Para visualizar os resultados no ANSYS MAPDL, use os seguintes comandos:")
        print(f"1. RESUME,'{model_name}_after_solve','db'")
        print("2. /POST1")
        print("3. SET,1,1")
        print("4. PLDISP")
        # Abrir a GUI manualmente
        print("Para abrir a GUI manualmente, execute o seguinte comando no terminal:")
        print(f"mapdl.open_gui()")
    
        # Encerrar o MAPDL corretamente
        mapdl.exit()