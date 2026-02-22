import numpy as np




def single_effect_model_derive():
    x_orig=np.array([1.0, 1.5, 2.0])
    x= x_orig - np.mean(x_orig)
    x_vec=x.reshape(-1,1)
    x_transpose = np.transpose(x_vec)
    sigma_0_squared=1.0
    sigma_squared=3.0
    y_orig=np.array(
        [
            2.1,
            2.9,
            3.8
        ]
    )
    y= y_orig - np.mean(y_orig)
    y_vec=y.reshape(-1,1)
    b_hat = np.dot(x,y)/np.dot(x,x)
    s_squared= sigma_squared/np.dot(x,x)
    sigma_1_squared =1/(1/s_squared + 1/sigma_0_squared)
    ratio = sigma_1_squared/s_squared
    mu_1 =sigma_1_squared/s_squared*b_hat
    derived_mean = sigma_0_squared*x_transpose @np.linalg.inv(sigma_0_squared*x_vec@x_transpose + sigma_squared*np.eye(3)   ) @ y_vec
    derived_mean_2 = x_transpose @np.linalg.inv(x_vec@x_transpose + sigma_squared/sigma_0_squared*np.eye(3)   ) @ y_vec
    derived_mean_3 = x_transpose@( sigma_0_squared/sigma_squared*np.eye(3)
                                   - (sigma_0_squared/sigma_squared)**2*(x_vec@x_transpose)/(1+(sigma_0_squared/sigma_squared)*x_transpose@x_vec))@ y_vec
    derived_mean_4 = sigma_0_squared/sigma_squared*(x_transpose@y_vec)/(1+sigma_0_squared/sigma_squared*(x_transpose@x_vec))
    derived_ratio = sigma_0_squared/sigma_squared / (1/np.dot(x,x,)+sigma_0_squared/sigma_squared)
    print(f"{b_hat=}, {s_squared=}, {sigma_1_squared=}, {mu_1=}, {derived_mean=},{derived_mean_2=}, {derived_mean_3=}, {derived_mean_4=}")
    print(f"{ratio=}, {derived_ratio=}")



if __name__ == '__main__':
    single_effect_model_derive()