import math
import numpy as np
from decimal import Decimal

#   input data
L = 39.2
w = 3
h = 0.18
d = 6.4
t_d = 0
v = 8.5
T_o = Decimal(292)
a = 2.48
Q = 410
eff = 0.3
q = Q*eff
k = 6.7E-03
n_layers = 63
x, y, z = 19.6, 0, -6.4
t_total = 300

#   start time of jth layer
t_start = [(((j-1) * L/v) + ((j-1) * t_d)) for j in range(n_layers)]

f = open('TC1_WALL1_M3.txt', 'w')

t = 1
while(t <= t_total):    
    for j in range(n_layers):
        if t_start[j] < t and t < t_start[j] + L/v:
            #   jth layer under deposition
            #   current layer = j
            delta_t = t - t_start[j]
            if j%2 != 0:
                real_x = v * (t - t_start[j])
                delta_x = x - real_x
            else:
                real_x = L - v * (t - t_start[j])
                delta_x = real_x - x
            
            real_z = (j - 0.5) * h
            delta_z = z - real_z
    
            #   original real heat source
            R_o = math.sqrt(delta_x**2 + y**2 + delta_z**2)
            A_o = Decimal(1 - math.erf((v * delta_t + R_o)/(2 * math.sqrt(a * delta_t))))
            B_o = Decimal(1 + math.erf((v * delta_t - R_o)/(2 * math.sqrt(a * delta_t))))
            TF_o = Decimal(0.5) * (A_o * np.exp(Decimal(v * R_o / a)) + B_o)
            real_T_o = Decimal((q / (2*math.pi*k*R_o))) * Decimal(math.exp(-v * (delta_x + R_o) / (2*a))) * TF_o
    
            #   bottom image of real heat source
            delta_z_bt = z + real_z + (2*d)
            R_bt = math.sqrt(delta_x**2 + y**2 + delta_z_bt**2)
            A_bt = Decimal(1 - math.erf((v * delta_t + R_bt)/(2 * math.sqrt(a * delta_t))))
            B_bt = Decimal(1 + math.erf((v * delta_t - R_bt)/(2 * math.sqrt(a * delta_t))))
            TF_bt = Decimal(0.5) * (A_bt * np.exp(Decimal(v * R_bt / a)) + B_bt)
            real_T_bt = Decimal((q / (2*math.pi*k*R_bt))) * Decimal(math.exp(-v * (delta_x + R_bt) / (2*a))) * TF_bt        

            #   total real heat source temp for jth layer
            real_T = real_T_o + real_T_bt
            
            #   calcs for i (= 1 to j-1) +ve and -ve virtual heat sources
            i = 1
            virtual_T = 0
            while(i < j):
                delta_t_pos = t - t_start[i]
                delta_t_neg = t - (t_start[i] + L/v)
                if i%2 != 0:
                    virtual_x = v * (t - t_start[i])
                    delta_x_virtual = x - virtual_x
                else:
                    virtual_x = L - v * (t - t_start[i])
                    delta_x_virtual = virtual_x - x
                virtual_z = (i - 0.5) * h
                delta_z_virtual = z - virtual_z
    
                #   original +ve virtual heat source
                R_o_pos = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_virtual**2)
                A_o_pos = Decimal(1 - math.erf((v * delta_t_pos + R_o_pos)/(2 * math.sqrt(a * delta_t_pos))))
                B_o_pos = Decimal(1 + math.erf((v * delta_t_pos - R_o_pos)/(2 * math.sqrt(a * delta_t_pos))))
                TF_o_pos = Decimal(0.5) * (A_o_pos * np.exp(Decimal(v * R_o_pos / a)) + B_o_pos)
                T_o_pos = Decimal((q / (2*math.pi*k*R_o_pos))) * Decimal(math.exp(-v * (delta_x_virtual + R_o_pos) / (2*a))) * TF_o_pos
        
                #   bottom image of +ve virtual heat source
                delta_z_bt_virtual = z + virtual_z + (2*d)
                R_bt_pos = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_bt_virtual**2)
                A_bt_pos = Decimal(1 - math.erf((v * delta_t_pos + R_bt_pos)/(2 * math.sqrt(a * delta_t_pos))))
                B_bt_pos = Decimal(1 + math.erf((v * delta_t_pos - R_bt_pos)/(2 * math.sqrt(a * delta_t_pos))))
                TF_bt_pos = Decimal(0.5) * (A_bt_pos * np.exp(Decimal(v * R_bt_pos / a)) + B_bt_pos)
                T_bt_pos = Decimal((q / (2*math.pi*k*R_bt_pos))) * Decimal(math.exp(-v * (delta_x_virtual + R_bt_pos) / (2*a))) * TF_bt_pos        
        
                #   total +ve virtual heat source temp for ith layer
                T_pos = T_o_pos + T_bt_pos
                    
                #   original -ve virtual heat source
                R_o_neg = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_virtual**2)
                A_o_neg = Decimal(1 - math.erf((v * delta_t_neg + R_o_neg)/(2 * math.sqrt(a * delta_t_neg))))
                B_o_neg = Decimal(1 + math.erf((v * delta_t_neg - R_o_neg)/(2 * math.sqrt(a * delta_t_neg))))
                TF_o_neg = Decimal(0.5) * (A_o_neg * np.exp(Decimal(v * R_o_neg / a)) + B_o_neg)
                T_o_neg = Decimal((-q / (2*math.pi*k*R_o_neg))) * Decimal(math.exp(-v * (delta_x_virtual + R_o_neg) / (2*a))) * TF_o_neg
        
                #   bottom image of -ve virtual heat source
                delta_z_bt_virtual = z + virtual_z + (2*d)
                R_bt_neg = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_bt_virtual**2)
                A_bt_neg = Decimal(1 - math.erf((v * delta_t_neg + R_bt_neg)/(2 * math.sqrt(a * delta_t_neg))))
                B_bt_neg = Decimal(1 + math.erf((v * delta_t_neg - R_bt_neg)/(2 * math.sqrt(a * delta_t_neg))))
                TF_bt_neg = Decimal(0.5) * (A_bt_neg * np.exp(Decimal(v * R_bt_neg / a)) + B_bt_neg)
                T_bt_neg = Decimal((-q / (2*math.pi*k*R_bt_neg))) * Decimal(math.exp(-v * (delta_x_virtual + R_bt_neg) / (2*a))) * TF_bt_neg        
        
                #   total -ve virtual heat source temp for ith layer
                T_neg = T_o_neg + T_bt_neg     
            
                virtual_T += T_pos + T_neg        
                i += 1
            
            total_T = T_o + real_T + virtual_T
            print(total_T - 273, file = f)
                 
    if t_start[n_layers - 1] + L/v < t:
        #   dwelling period
        #   previous layer = n_layers - 1    
            
        #   calcs for i (1 to n_layers - 1) +ve and -ve virtual heat sources
        i = 1
        virtual_T = 0
        while(i <= n_layers - 1):
            delta_t_pos = t - t_start[i]
            delta_t_neg = t - (t_start[i] + L/v)
            if i%2 != 0:
                virtual_x = v * (t - t_start[i])
                delta_x_virtual = x - virtual_x
            else:
                virtual_x = L - v * (t - t_start[i])
                delta_x_virtual = virtual_x - x
            
            virtual_z = (i - 0.5) * h
            delta_z_virtual = z - virtual_z
    
            #   original +ve virtual heat source
            R_o_pos = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_virtual**2)
            A_o_pos = Decimal(1 - math.erf((v * delta_t_pos + R_o_pos)/(2 * math.sqrt(a * delta_t_pos))))
            B_o_pos = Decimal(1 + math.erf((v * delta_t_pos - R_o_pos)/(2 * math.sqrt(a * delta_t_pos))))
            TF_o_pos = Decimal(0.5) * (A_o_pos * np.exp(Decimal(v * R_o_pos / a)) + B_o_pos)
            T_o_pos = Decimal((q / (2*math.pi*k*R_o_pos))) * Decimal(math.exp(-v * (delta_x_virtual + R_o_pos) / (2*a))) * TF_o_pos
    
            #   bottom image of +ve virtual heat source
            delta_z_bt_virtual = z + virtual_z + (2*d)
            R_bt_pos = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_bt_virtual**2)
            A_bt_pos = Decimal(1 - math.erf((v * delta_t_pos + R_bt_pos)/(2 * math.sqrt(a * delta_t_pos))))
            B_bt_pos = Decimal(1 + math.erf((v * delta_t_pos - R_bt_pos)/(2 * math.sqrt(a * delta_t_pos))))
            TF_bt_pos = Decimal(0.5) * (A_bt_pos * np.exp(Decimal(v * R_bt_pos / a)) + B_bt_pos)
            T_bt_pos = Decimal((q / (2*math.pi*k*R_bt_pos))) * Decimal(math.exp(-v * (delta_x_virtual + R_bt_pos) / (2*a))) * TF_bt_pos        
    
            #   total +ve virtual heat source temp for ith layer
            T_pos = T_o_pos + T_bt_pos
                
            #   original -ve virtual heat source
            R_o_neg = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_virtual**2)
            A_o_neg = Decimal(1 - math.erf((v * delta_t_neg + R_o_neg)/(2 * math.sqrt(a * delta_t_neg))))
            B_o_neg = Decimal(1 + math.erf((v * delta_t_neg - R_o_neg)/(2 * math.sqrt(a * delta_t_neg))))
            TF_o_neg = Decimal(0.5) * (A_o_neg * np.exp(Decimal(v * R_o_neg / a)) + B_o_neg)
            T_o_neg = Decimal((-q / (2*math.pi*k*R_o_neg))) * Decimal(math.exp(-v * (delta_x_virtual + R_o_neg) / (2*a))) * TF_o_neg
    
            #   bottom image of -ve virtual heat source
            delta_z_bt_virtual = z + virtual_z + (2*d)
            R_bt_neg = math.sqrt(delta_x_virtual**2 + y**2 + delta_z_bt_virtual**2)
            A_bt_neg = Decimal(1 - math.erf((v * delta_t_neg + R_bt_neg)/(2 * math.sqrt(a * delta_t_neg))))
            B_bt_neg = Decimal(1 + math.erf((v * delta_t_neg - R_bt_neg)/(2 * math.sqrt(a * delta_t_neg))))
            TF_bt_neg = Decimal(0.5) * (A_bt_neg * np.exp(Decimal(v * R_bt_neg / a)) + B_bt_neg)
            T_bt_neg = Decimal((-q / (2*math.pi*k*R_bt_neg))) * Decimal(math.exp(-v * (delta_x_virtual + R_bt_neg) / (2*a))) * TF_bt_neg        
    
            #   total -ve virtual heat source temp for ith layer
            T_neg = T_o_neg + T_bt_neg
            
            virtual_T += T_pos + T_neg        
            i += 1
            
        total_T = T_o + virtual_T
        print(total_T - 273, file = f)       
    t += 1
    
f.close()