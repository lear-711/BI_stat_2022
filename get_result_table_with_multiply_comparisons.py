import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests


def get_result_table(first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table, multiply_comp_method = False, adj_alpha = 0.05):

    # Находим таблицы по путям, считываем их и записываем в переменную
    first_data_path = first_cell_type_expressions_path
    second_data_path = second_cell_type_expressions_path
    first_table = pd.read_csv(first_data_path, index_col=0)
    second_table = pd.read_csv(second_data_path, index_col=0)


    # Считаем, значимо ли отличается средняя экспрессия генов между клеточными типами:
    def check_dge_with_ci(first_table, second_table):
    # dge - differential gene expression
        ci_test_results = []
        for gene in first_table:
          if (first_table[gene].dtype == 'float64') or (first_table[gene].dtype == 'int'):
            interval_1 = st.t.interval(alpha=0.95,
                                      df=len(first_table[gene]) - 1,
                                      loc=np.mean(first_table[gene]),
                                      scale=st.sem(first_table[gene]))
            interval_2 = st.t.interval(alpha=0.95,
                                      df=len(second_table[gene]) - 1,
                                      loc=np.mean(second_table[gene]),
                                      scale=st.sem(second_table[gene]))
            ci_test_results.append(check_intervals_intersect(interval_1, interval_2))
        return ci_test_results
      
    ci_test_results_d = check_dge_with_ci(b_cells_expression_data, nk_cells_expression_data)


    # Используем z-критерий для расчета отличий между средними экспрессиями генов между клеточными линиями:
    def check_dge_with_ztest(first_table, second_table):
    # dge - differential gene expression
        z_test_results = []
        z_test_p_values = []

        for gene in first_table:
            if (first_table[gene].dtype == 'float64') or (first_table[gene].dtype == 'int'):
                z_test = ztest(
                    first_table[gene],
                    second_table[gene]
                    )
                z_test_p_values.append(z_test[1])
                if z_test[1] < 0.05:
                    z_test_results.append(False)
                else:
                    z_test_results.append(True)
        return z_test_results, z_test_p_values

    z_test_results_d, z_test_p_values_d = check_dge_with_ztest(b_cells_expression_data, nk_cells_expression_data)


    # Если передано значение для аргумента multiply_comp_method, отличное от False, считаем с поправкой на множественные сравнения:
    if multiply_comp_method != False:
      corr_p_values = multipletests(z_test_p_values_d, adj_alpha, method = multiply_comp_method)


    # Находим разницу в средних экспрессиях между 1 и 2 таблицами для каждого гена:
    def calc_mean_difference(first_table, second_table):

        mean_diff = []

        for gene in first_table:
          if (first_table[gene].dtype == 'float64') or (first_table[gene].dtype == 'int'):
            mean_d = np.mean(first_table[gene]) - np.mean(second_table[gene])
            mean_d = round(mean_d, 2)
            mean_diff.append(mean_d)
        return mean_diff

    mean_diff_results = calc_mean_difference(b_cells_expression_data, nk_cells_expression_data)


    # Сохраняем в список названия генов:
    gene_names = []
    for gene in expression_data:
      gene_names.append(gene)
    gene_names.pop()


    # Создаем словарь с результатами
    if multiply_comp_method == False: 
      results = {
      "gene_names": gene_names,
      "ci_test_results": ci_test_results_d,
      "z_test_p_values": z_test_p_values_d,
      "z_test_results": z_test_results_d,
      "mean_diff": mean_diff_results
      }
    else:
      results = {
      "gene_names": gene_names,
      "z_test_p_values": z_test_p_values_d,
      "multiply_adj_results": corr_p_values[0],
      "adjusted_p_values": corr_p_values[1],
      "mean_diff": mean_diff_results
      }
    for k in results:
      print(k, len(results[k]))
    # Из словаря делаем датафрейм
    results = pd.DataFrame(results)
    results.head()

    results.to_csv("/content/drive/MyDrive/Статистика ИБ/" + save_results_table)