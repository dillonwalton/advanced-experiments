import math
import statistics

def model(
        x,

        rb87_f2_a,
        rb87_f2_x0,
        rb87_f2_sigma,

        rb85_f3_a,
        rb85_f3_x0,
        rb85_f3_sigma,

        rb85_f2_a,
        rb85_f2_x0,
        rb85_f2_sigma,

        rb87_f1_a,
        rb87_f1_x0,
        rb87_f1_sigma,

        lin_m,
        lin_b
):
    gaussian_1 = rb87_f2_a*math.exp(-(x-rb87_f2_x0)**2/(2*rb87_f2_sigma**2))
    gaussian_2 = rb85_f3_a*math.exp(-(x-rb85_f3_x0)**2/(2*rb85_f3_sigma**2))
    gaussian_3 = rb85_f2_a*math.exp(-(x-rb85_f2_x0)**2/(2*rb85_f2_sigma**2))
    gaussian_4 = rb87_f1_a*math.exp(-(x-rb87_f1_x0)**2/(2*rb87_f1_sigma**2))
    bg = lin_m * x + lin_b
    return gaussian_1 + gaussian_2 + gaussian_3 + gaussian_4 + bg

def temperature(sigma, f_0, isotope):

    c = 3 * (10**8)
    k_B = 1.3807 * (10**(-23))
    if isotope == "85Rb":
        M = 1.409993 * (10**(-25)) # https://steck.us/alkalidata/rubidium85numbers.pdf
    elif isotope == "87Rb":
        M = 1.443160 * (10**(-25)) # https://steck.us/alkalidata/rubidium87numbers.pdf
    
    return ((sigma * c)/f_0)**2 * (M / k_B)

def uncertainty(sigma, sigma_sigma, f_0, sigma_f_0, isotope):

    c = 3 * (10**8)
    k_B = 1.3807 * (10**(-23))
    if isotope == "85Rb":
        M = 1.409993 * (10**(-25)) # https://steck.us/alkalidata/rubidium85numbers.pdf
    elif isotope == "87Rb":
        M = 1.443160 * (10**(-25)) # https://steck.us/alkalidata/rubidium87numbers.pdf

    f_0_weight = ((2 * c**2 * sigma**2 * M)/(k_B * f_0**3))**2
    sigma_weight = ((2 * sigma * c**2 * M)/(k_B * f_0**2))**2

    return math.sqrt((f_0_weight * sigma_f_0**2) + (sigma_weight * sigma_sigma**2))

def weighted_avg(data):
    # weight
    weight_sum = 0
    j = 0
    for j in range(0, len(data)):
        weight = (1/((data[j][1])**2))
        weight_sum += weight
        j += 1
    
    # average
    average_sum = 0
    i = 0
    for i in range(0, len(data)):
        normalized_weight = (1/((data[i][1])**2)) / weight_sum
        value = data[i][0]

        average_sum += (value * normalized_weight)

        i+=1
    
    # uncert
    temp_std = statistics.stdev([data[0][0],data[1][0],data[2][0],data[3][0]])
    weighted_variance_sum = 0
    y = 0
    for y in range(0, len(data)):
        normalized_weight = (1/((data[y][1])**2)) / weight_sum
        weighted_variance_sum += (normalized_weight**2)
    weighted_avg_uncert = temp_std * math.sqrt(weighted_variance_sum)
    
        
    return average_sum,weighted_avg_uncert

