

I = np.arange(10)

mask = np.array([1,1,1,1,0,0,0,0,0,0]).astype('bool')

J = I.copy()
#np.random.seed(0)
J[mask] = np.random.permutation(J[mask])
print J
