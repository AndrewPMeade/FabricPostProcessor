import math


class Distribution():
	def __init__(self, pvals):
		self.Paramters = pvals
		self.Name = None

	def __str__(self):
		return f"{self.Name} {*self.Paramters,}"

	def PDF(self, x):
		pass

	def CDF(self, x):
		pass


class Weibull(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "Weibull"
		assert(len(pvals) == 2)

	def PDF(self, X):
		L = self.Paramters[0]
		K = self.Paramters[1]

		Ret = -math.pow(X/L, K)
		Ret = math.pow(math.e, Ret)

		Ret = (K/L) * math.pow(X/L, K-1) * Ret

		return Ret

	def CDF(self, X):
		L = self.Paramters[0]
		K = self.Paramters[1]

		return 1.0 - math.pow(math.e, -math.pow(X/L, K))

class LogNormal(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "LogNormal"
		assert(len(pvals) == 2)

	def PDF(self, X):
		M = self.Paramters[0]
		S = self.Paramters[1]

		Ret = (math.log(X) - M)**2 / (2.0 * S**2)
		Ret = math.exp(-Ret)
		return (1.0 / (X * S * math.sqrt(2.0 * math.pi))) * Ret


	def CDF(self, X):
		M = self.Paramters[0]
		S = self.Paramters[1]

		Ret = 1.0  + math.erf((math.log(X) - M) / (S * math.sqrt(2.0)))

		return 0.5 * Ret

class Gamma(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "Gamma"
		assert(len(pvals) == 2)

	def PDF(self, X):
		Kappa = self.Paramters[0]
		Theta = self.Paramters[1]

		A = 1.0 / (math.gamma(Kappa) * math.pow(Theta, Kappa))
		B = math.pow(X, Kappa-1.0)
		C = math.exp(-(X/Theta))

		return A * B * C

	def CDF(self, X):
		pass

class Uniform(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "Uniform"
		assert(len(pvals) == 2)

	def PDF(self, X):
		Min = self.Paramters[0]
		Max = self.Paramters[1]

		if X < Min or X > Max:
			return 0.0;

		return 1.0 / (Max - Min)

	def CDF(self, X):
		Min = self.Paramters[0]
		Max = self.Paramters[1]

		if X < Min:
			return 0.0

		if X > Max:
			return 1.0

		return (X - Min) / (Max - Min)

class Chi(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "Chi"
		assert(len(pvals) == 1)

	def PDF(self, X):
		Kappa = self.Paramters[0]
		
		A = 1.0 / (math.pow(2.0, Kappa / 2.0) * math.gamma(Kappa / 2.0))
		B = math.pow(X, (Kappa / 2.0) -1.0)
		C = math.exp(-X/2)

		return A * B *C 

	def CDF(self, X):
		pass

# NB the beta value has not effect, possibly due to scaling nature. 
class SGamma(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "SGamma"
		assert(len(pvals) == 2)

	def PDF(self, X):
		Alpha = self.Paramters[0]
		Beta = self.Paramters[1]

		S = 1.0 / ((Alpha - 1.0) * Beta)
		A = math.exp(-(X/S) / Beta)
		B = math.pow(X/S, -1 + Alpha)
		B = A * B * math.pow(Beta, -Alpha)

		Ret = (B / math.gamma(Alpha)) / S

		return Ret

	
	def CDF(self, X):
		pass

		
class Exp(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)
		
		self.Name = "Exp"
		assert(len(pvals) == 1)

	def PDF(self, X):
		Lambda = 1.0 / self.Paramters[0]		

		return Lambda * math.exp(-Lambda * X)

	def CDF(self, X):
		pass

# normal, paramter is sig not sig2
class Normal(Distribution):
	def __init__(self, pvals):
		super().__init__(pvals)

		self.Name = "Normal"
		assert(len(pvals) == 2)

	def PDF(self, X):

		Mu = self.Paramters[0]
		Sig = self.Paramters[1]
		
		A = 1.0 / (Sig * math.sqrt(2.0 * math.pi))
		B = ((X - Mu) / Sig) ** 2
		B = math.exp(-0.5 * B)

		return A * B

	
	def CDF(self, X):
		pass



def CreateDistribution(Name, pvals):

	Name = Name.lower()

	if Name == "weibull":
		return Weibull(pvals)

	if Name == "lognormal":
		return LogNormal(pvals)

	if Name == "gamma":
		return Gamma(pvals)

	if Name == "uniform":
		return Uniform(pvals)

	if Name == "chi":
		return Chi(pvals)

	if Name == "exp":
		return Exp(pvals)

	if Name == "sgamma":
		return SGamma(pvals)

	if Name == "normal":
		return Normal(pvals)

	return None