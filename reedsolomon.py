import copy


class ReedSolomon:
	# Galois fields
	# -- exponents (anti-logarithms)
	__GFEXP = [0] * 512
	
	# -- logarithms
	__GFLOG = [0] * 256
	
	
	def __init__(self):
		# prepare the exponential and logarithmic fields
		self.__GFEXP[0] = 1
		byteValue = 1
		for bytePos in range(1, 255):
			byteValue <<= 1
			if (byteValue & 0x100):
				byteValue ^= 0x11d
			
			# update the field elements
			self.__GFEXP[bytePos] = byteValue
			self.__GFLOG[byteValue] = bytePos
		
		# finalise the exponential field
		for bytePos in range(255,512):
			self.__GFEXP[bytePos] = self.__GFEXP[bytePos - 255]
	
	
	
	## GALOIS PRIMITIVE OPERATIONS
	# -----
	# Galois multiplication
	# argX, argY: multiplicand, multiplier
	# byteValue: product
	def __gfMult(self, argX, argY):
		# parametre checks
		if ((argX == 0) or (argY == 0)):
			byteValue = 0
		else:
			# perform the operation
			byteValue = self.__GFLOG[argX]
			byteValue += self.__GFLOG[argY]
			byteValue = self.__GFEXP[byteValue]
		
		# return the product result
		return (byteValue)
	
	# Galois division
	# argX, argY: dividend, divisor
	# byteValue: quotient
	def __gfDivi(self, argX, argY):
		# validate the divisor
		if (argY == 0):
			raise ZeroDivisionError()
			
		# validate the dividend
		if (argX == 0):
			byteValue = 0
		else:
			# perform the division
			byteValue = self.__GFLOG[argX] - self.__GFLOG[argY]
			byteValue += 255
			byteValue = self.__GFEXP[byteValue]
		
		# return the division result
		return (byteValue)
	
	
	## GALOIS POLYNOMIAL OPERATIONS
	# -----
	# Polynomial addition
	# polyA, polyB: polynomial addends
	# polySum: polynomial sum
	def _gfPolyAdd(self, polyA, polyB):
		# initialise the polynomial sum
		polySum = [0] * max(len(polyA), len(polyB))
		
		# process the first addend
		for polyPos in range(0, len(polyA)):
			polySum[polyPos + len(polySum) - len(polyA)] = polyA[polyPos]
		
		# add the second addend
		for polyPos in range(0, len(polyB)):
			polySum[polyPos + len(polySum) - len(polyB)] ^= polyB[polyPos]
		
		# return the sum
		return (polySum)
		
	
	# Polynomial multiplication
	# polyA, polyB: polynomial factors
	# polyProd: polynomial product
	def _gfPolyMult(self, polyA, polyB):
		# initialise the product
		polyProd = len(polyA) + len(polyB) - 1
		polyProd = [0] * polyProd
		
		# start multiplying
		for posB in range(0, len(polyB)):
			for posA in range(0, len(polyA)):
				polyProd[posA + posB] ^= self.__gfMult(polyA[posA], polyB[posB])
		
		# return the product result
		return (polyProd)
	
	
	# Polynomial scaling
	# argPoly: polynomial argument
	# argX: scaling factor
	# polyVal: scaled polynomial
	def _gfPolyScale(self, argPoly, argX):
		# initialise the scaled polynomial
		polyVal = [0] * len(argPoly)
		
		# start scaling
		for polyPos in range(0, len(argPoly)):
			polyVal[polyPos] = self.__gfMult(argPoly[polyPos], argX)
		
		# return the scaled polynomial
		return (polyVal)
	
	
	# Polynomial eValueation
	# argPoly: polynomial argument
	# argX: independent variable
	# byteValue: dependent variable
	def _gfPolyEval(self, argPoly, argX):
		# initialise the polynomial result
		byteValue = argPoly[0]
		
		# eValueate the polynomial argument
		for polyPos in range(1, len(argPoly)):
			tempValue = self.__gfMult(byteValue, argX) 
			tempValue = tempValue ^ argPoly[polyPos]
			byteValue = tempValue
		
		# return the eValueated result
		return (byteValue)
	
	
	
	## REED-SOLOMON SUPPORT ROUTINES
	# -----
	# Prepare the generator polynomial
	# errSize: number of error symbols
	# polyValue: generator polynomial
	def _rsGenPoly(self, errSize):
		polyValue = [1]
		
		for polyPos in range(0, errSize):
			tempVal = [1, self.__GFEXP[polyPos]]
			polyValue = self._gfPolyMult(polyValue, tempVal)
		
		# return the polynomial result
		return (polyValue)
	
	
	
	## REED-SOLOMON ENCODING
	# ------
	# argMesg: the message block
	# errSize: number of error symbols
	# outBuffer: the encoded output buffer
	def RSEncode(self, argMesg, errSize):
		
		# prepare the generator polynomial
		polyGen = self._rsGenPoly(errSize)
		
		# prepare the output buffer
		outBuffer = (len(argMesg) + errSize)
		outBuffer = [0] * outBuffer
		
		# initialise the output buffer
		for mesgPos in range(0, len(argMesg)):
			mesgChar = argMesg[mesgPos]
			outBuffer[mesgPos] = ord(mesgChar)
		
		# begin encoding
		for mesgPos in range(0, len(argMesg)):
			mesgChar = outBuffer[mesgPos]
			if (mesgChar != 0):
				for polyPos in range(0, len(polyGen)):
					tempValue = self.__gfMult(polyGen[polyPos], mesgChar)
					outBuffer[mesgPos + polyPos] ^= tempValue
		
		# finalise the output buffer
		for mesgPos in range(0, len(argMesg)):
			mesgChar = argMesg[mesgPos]
			outBuffer[mesgPos] = ord(mesgChar)
		
		# return the output buffer
		return (outBuffer)
	
	
	## REED-SOLOMON DECODING
	# -----
	# Generate the syndrome polynomial
	# argCode: the code block
	# errSize: number of error symbols
	# polyValue: the syndrome polynomial
	def _rsSyndPoly(self, argCode, errSize):
		# initialise the polynomial
		polyValue = [0] * errSize
		
		# compute the polynomial terms
		for errPos in range(0, errSize):
			byteValue = self.__GFEXP[errPos] 
			polyValue[errPos] = self._gfPolyEval(argCode, byteValue)
		
		# return the polynomial
		return (polyValue)
	
	# The Forney algorithm
	# polySynd: the syndrome polynomial
	# eraseLoci: list of erasures
	# errSize: number of error symbols
	# polyValue: the error locator polynomial 
	def _rsForney(self, polySynd, eraseLoci, errSize):
		# make a copy of the syndrome polynomial
		polyValue = list(polySynd)
		
		# compute the polynomial terms
		for posI in range(0, len(eraseLoci)):
			termX = errSize - 1 - eraseLoci[posI]
			termX = self.__GFEXP[termX]
			for posJ in range(0, len(polyValue) - 1):
				termY = self.__gfMult(polyValue[posJ], termX)
				termY ^= polyValue[posJ + 1]
				polyValue[posJ] = termY
			polyValue.pop()
		
		# return the polynomial result
		return (polyValue)
	
	# Locate the message errors
	# errLoci: error locator polynomial
	# errSize: number of error symbols
	def _rsFindErr(self, errLoci, errSize):
		# initialise the polynomial locals
		errPoly = [1]
		tempPoly = [1]
		
		# generate the error locator polynomial
		# - Berklekamp-Massey algorithm
		for posSynd in range(0, len(errLoci)):
			tempPoly.append(0)
			termSynd = errLoci[posSynd]
			
			for posErr in range(1, len(errPoly)):
				termPoly = errPoly[len(errPoly) - posErr - 1]
				termPoly = self.__gfMult(termPoly, errLoci[posSynd - posErr])
				termSynd ^= termPoly
			
			if (termSynd != 0):
				if (len(tempPoly) > len(errPoly)):
					tNewP = self._gfPolyScale(tempPoly, termSynd)
					tempPoly = self._gfPolyScale(errPoly, self.__gfDivi(1, termSynd))
					errPoly = tNewP
				
				tempValue = self._gfPolyScale(tempPoly, termSynd)
				errPoly = self._gfPolyAdd(errPoly, tempValue)
		
		# count the number of errors
		errCount = len(errPoly) - 1
		if ((errCount * 2) > len(errLoci)):
			print("Zbyt wiele błędów do poprawienia")
			return (None)
		else:
			print("Liczba błędów: ", errCount)
		
		# calculate the polynomial zeroes
		errList = []
		for errPos in range(0, errSize):
			errZed = self._gfPolyEval(errPoly, self.__GFEXP[255 - errPos])
			if (errZed == 0):
				errZed = errSize - errPos - 1
				errList.append(errZed)
		
		if (len(errList) != errCount):
			print("Nie udało się zlokalizować błędów")
			return (None)
		else:
			return (errList)
	
	# Correct errors and erasures
	# argCode: the message code block
	# polySynd: the sydrome polynomial
	# errList: list of error and erasure positions
	def _rsCorrect(self, argCode, polySynd, errList):
		# prepare the locator polynomial
		polyLoci = [1]
		for errPos in range(0, len(errList)):
			errTerm = len(argCode) - errList[errPos] - 1
			errTerm = self.__GFEXP[errTerm]
			polyLoci = self._gfPolyMult(polyLoci, [errTerm, 1])
		
		# prepare the error eValueator polynomial
		errEval = polySynd[0:len(errList)]
		errEval.reverse()
		errEval = self._gfPolyMult(errEval, polyLoci)
		
		tMark = len(errEval) - len(errList)
		errEval = errEval[tMark:len(errEval)]
		
		# the error locator polynomial, minus even terms
		errLoci = polyLoci[len(polyLoci) % 1 : len(polyLoci) : 2]
		
		# start correcting
		for errPos in range(0, len(errList)):
			errByte = errList[errPos] - len(argCode) + 256
			errByte = self.__GFEXP[errByte]
			
			errValue = self._gfPolyEval(errEval, errByte)
			
			errAdj = self.__gfMult(errByte, errByte)
			errAdj = self._gfPolyEval(errLoci, errAdj)
			
			mesgByte = self.__gfMult(errByte, errAdj)
			mesgByte = self.__gfDivi(errValue, mesgByte)
			argCode[errList[errPos]] ^= mesgByte
		return (argCode)
	
	# Main decode routine
	# argCode: the message code block
	# errSize: number of error symbols
	def RSDecode(self, argCode, errSize):
		
		# initialise the code buffer
		codeBuffer = list(argCode)
		
		# count the number of erasures
		eraseCount = []
		for codePos in range(0, len(codeBuffer)):
			if (codeBuffer[codePos] < 0):
				codeBuffer[codePos] = 0
				eraseCount.append(codePos)
		if (len(eraseCount) > errSize):
			print("Zbyt wiele wymazań")
			return (None)
		
		# prepare the syndrome polynomial
		polySynd = self._rsSyndPoly(codeBuffer, errSize)
		if (max(polySynd) == 0):
			print("Wiadomość nie zawiera błędów")
			return (codeBuffer)
		
		# prepare the error locator polynomial
		errLoci = self._rsForney(polySynd, eraseCount, len(codeBuffer))
		
		# locate the message errors
		errList = self._rsFindErr(errLoci, len(codeBuffer))
		if (errList == None):
			print("Nie znaleziono żadnych błędów")
			return (None)
		else:
			print("Zlokalizowane błędy: ", errList)
		
		# start correcting errors and erasures
		outMesg = self._rsCorrect(codeBuffer, polySynd, (eraseCount + errList))
		return (outMesg)



if __name__ == '__main__':
	# ------
	# MAIN TEST SCRIPT
	# ------
	fooRS = ReedSolomon()

	# set the test parametres
	tMesg = "Hello world!"
	# number of the redundancy symbols
	tSize = 6

	# encode the message
	tCode = fooRS.RSEncode(tMesg, tSize)
	print("Słowo kodowe: ", tCode)
	print("\r\r")

	# introduce errors/erasures
	errCode = copy.deepcopy(tCode)
	errCode[3] = 9
	errCode[7] = -1
	errCode[10] = 50
	print("Słowo kodowe (z 3 błędami/wymazaniami): ", errCode)
	print("\r\r")

	# decode the message
	decodedMesg = fooRS.RSDecode(errCode, tSize)
	print("Zdekodowana wiadomość: ", decodedMesg)
	try:
		print("Zdekodowany tekst: ", ''.join([chr(char) for char in decodedMesg][:len(decodedMesg) - tSize]))
	except TypeError:
		pass

	# print(tCode)
	# print(decodedMesg)
			
	assert tCode == decodedMesg
