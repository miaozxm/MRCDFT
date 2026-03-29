import numpy as np
import pandas as pd
import os

def HWG_Solution_f(HH, NN,option,cut,if_print_energy=False):
    """
    Print out the wave function and energy of the first state.  
    -HH: the Hamiltonian kernel matrix of deformation configurations.
    -NN: the Norm kernel matrix of deformation configurations.
    -option: choose 'cutoff' or 'nos'.
    -cut: cutoff parameter, commonly use 10**(-4) or 10**(-5).
    -if_print_energy: choose 'yes' or 'no'.
    """
    vals, vecs = np.linalg.eigh(NN)
    if option=='cutoff':
        basis = vecs[:, vals/vals[-1] >= cut]
        scale = vals[vals/vals[-1] >= cut]
        print("number of nature state:", len(scale))
        new_basis = basis / np.sqrt(scale)[np.newaxis, :]
        matrix_natural = np.linalg.multi_dot((new_basis.T, HH, new_basis))
        E, gg = np.linalg.eigh(matrix_natural)
        
        Energy=E[0]
        ff=np.dot(new_basis,gg)
        wfsff=ff[:,0]
        if if_print_energy=='yes':
            print("The energy of the ground state:", Energy)

    elif option == 'nos':
        basis = vecs[:,-cut:]
        scale = vals[-cut:]
        print("number of nature state:", len(scale))
        new_basis = basis / np.sqrt(scale)[np.newaxis, :]
        matrix_natural = np.linalg.multi_dot((new_basis.T, HH, new_basis))
        E, gg = np.linalg.eigh(matrix_natural)
        Energy=E[0]
        ff=np.dot(new_basis,gg)
        wfsff=ff[:,0]
        if if_print_energy=='yes':
            print("The energy of the ground state:", Energy)
    # return Energy,wfsff
    return E, ff



def HWG_Solution_G(HH, NN,option,cut,if_print_energy=False):
    """
    Convert the mixing coefficients f to g.
    The input para are the same as Function HWG_Solution_f.
    """
    vals, vecs = np.linalg.eigh(NN)
    if option=='cutoff':
        basis = vecs[:, vals/vals[-1] >= cut]
        scale = vals[vals/vals[-1] >= cut]
        new_basis = basis / np.sqrt(scale)[np.newaxis, :]
        matrix_natural = np.linalg.multi_dot((new_basis.T, HH, new_basis))
        E, gg = np.linalg.eigh(matrix_natural)
        GG=np.dot(basis,gg)
        wfsGG=GG[:,0]
    elif option == 'nos':
        basis = vecs[:,-cut:]
        scale = vals[:-cut:]
        new_basis = basis / np.sqrt(scale)[np.newaxis, :]
        matrix_natural = np.linalg.multi_dot((new_basis.T, HH, new_basis))
        E, gg = np.linalg.eigh(matrix_natural)
        Energy=E[0]
        GG=np.dot(basis,gg)
        wfsGG=GG[:,0]
        if if_print_energy=='yes':
            print("The energy of the state:", Energy)
    # return wfsGG
    return E, GG


def add_kernel(filename,line_number,start_pos,end_pos,numbers_list):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            if 1 <= line_number <= len(lines):
                line = lines[line_number - 1]
                if 1 <= start_pos <= end_pos <= len(line):
                    substring = line[start_pos - 1:end_pos]
                    number = float(substring) 
                    numbers_list.append(number)
                else:
                    print("The start point or end point is beyond the range of line.")
            else:
                print("The line number is beyond the range of file.")
    except FileNotFoundError:
        print(f"Can not find the File {filename}.")
    except ValueError:
        print("The extracted string cannot be converted to a number.")
    except Exception as e:
        print(f"Error: {e}")


def generate_kernel(beta2s,kernel_path,kernel_name,line_number,if_print):
    """
    This function is for generating the diagonal matrix elements.
    Use the function generate_kernel, you can obtain 2 lists: numberlist_norm, numberlist_hamtotal.
    In numberlist_norm, the first dimension is different angular momentum, the second is different configurations, 
    the first and second number in the third dimension of the list are the beta2 of initial state and final state, 
    and the rest 1 number in the third dimension of the list is the Norm kernel.
    -beta2s: ["+00", "+01", "+02", "+03", "+04", "+05", "+06", "-01", "-02", "-03"] # change for different deformation configuration sets.
    -kernel_num: 55, matches (lens of beta2s)*(lens of beta2s + 1)/2
    -if_print: True or False.
    """
    kernel_num = int(len(beta2s) * (len(beta2s)+1)/2) 
    numberlist_norm = [[] for _ in range(kernel_num)]
    numberlist_hamtotal = [[] for _ in range(kernel_num)]
    # the number 5 matches the lens of line_number.
    # this is for generate the total hamiltonian kernel, q0 kernel and be2 kernel
    beta2i = 0
    k = 0
    for beta21 in beta2s:
        beta2f = 0
        for beta22 in beta2s:
            if beta2f >= beta2i:
                filename = os.path.join(kernel_path,kernel_name.format(beta21,beta22))
                # add the Hamiltonian kernel into the list numberlist_hamtotal
                # add_kernel(filename,2,1,15,numberlist_hamtotal[j][k])
                # add_kernel(filename,2,31,45,numberlist_hamtotal[j][k])
                numberlist_hamtotal[k].append(float(beta21)/10)
                numberlist_hamtotal[k].append(float(beta22)/10)
                add_kernel(filename,line_number,31,45,numberlist_hamtotal[k])
                # add the Norm kernel into the list numberlist_norm
                # add_kernel(filename,2,1,15,numberlist_norm[j][k])
                # add_kernel(filename,2,31,45,numberlist_norm[j][k])
                numberlist_norm[k].append(float(beta21)/10)
                numberlist_norm[k].append(float(beta22)/10)
                add_kernel(filename,line_number,1,15,numberlist_norm[k])
                numberlist_hamtotal[k][2] = numberlist_hamtotal[k][2]*numberlist_norm[k][2]
                k += 1
            beta2f += 1
        beta2i += 1

    if if_print:
        print(numberlist_norm)
        print(numberlist_hamtotal)
    return(numberlist_norm,numberlist_hamtotal)
    


def reshape_matrix(data_list,if_print, valuecols=0):
    """
    This function transforms a list of data into a symmetrical DataFrame and then extracts the numeric data as a 2D array. 
    The data list contains entries with row and column indices, 
    along with associated values and additional irrelevant information.
    -data_list: numberlist_totalham[J] or numberlist_norm[J], or numberlist_q0(be2)[J].
    -if_print: True or False.
    -valuecols: 1,2 3, ...
    """
    # Extracting the required row labels, column labels, and values
    rows = [row[0] for row in data_list]
    cols = [col[1] for col in data_list]
    values = [value[valuecols+2] for value in data_list]
    # Creating DataFrame
    df = pd.DataFrame({'row': rows, 'column': cols, 'value': values})
    # Using pivot_table to create a DataFrame with 'row' and 'column' as index and column labels
    df_pivot = df.pivot_table(index='row', columns='column', values='value')
    # Aligning the DataFrame's row and column indexes, ensuring the matrix is symmetric
    df_pivot = df_pivot.reindex(index=sorted(df_pivot.index.union(df_pivot.columns)),
                                columns=sorted(df_pivot.columns.union(df_pivot.index)))
    # Filling missing values at symmetric positions
    df_filled = df_pivot.combine_first(df_pivot.T)
    # Extracting the values matrix as a 2D array
    values_matrix = df_filled.values
    if if_print:
        print(values_matrix)
    return(values_matrix)
