import numpy as np 


x-col = np.array(faraz.ColA)

b = np.array(faraz['ColA'])



from sklean.linear_model import LinearRegression 


x_train, y_train = faraz[['ColA']], faraz[['ColB']]

faraz_regr = LinearRegression()
faraz_regr.fit(x_train, y_train)