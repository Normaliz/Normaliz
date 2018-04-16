# encoding=utf8

import PyNormaliz_cpp
from PyNormaliz_cpp import *


class Cone:

    def __init__(self, *args, **kwargs):
        input_list = [k for k in args]
        for i in kwargs:
            current_input = kwargs[i]
            if type(current_input) == list and len(current_input) > 0 and type(current_input[0]) != list:
                kwargs[i] = [current_input]
            elif type(current_input) == bool and current_input == True:
                kwargs[i] = current_input = [[]]
            elif type(current_input) == bool and current_input == False:
                kwargs.pop(i)
        self.cone = PyNormaliz_cpp.NmzCone(input_list,**kwargs)

    def __process_keyword_args(self, keywords):
        input_list = []
        for i in keywords:
            if keywords[i] == True:
                input_list.append(i)
        return input_list

    def print_properties(self):
        props = PyNormaliz_cpp.NmzListConeProperties()
        goals = props[0]
        for x in goals:
            if (PyNormaliz_cpp.NmzIsComputed(self.cone, x)):
                print(x + ":")
                print(PyNormaliz_cpp.NmzResult(self.cone, x))
                print("\n")

    def __str__(self):
        return "<Normaliz Cone>"

    def __repr__(self):
        return "<Normaliz Cone>"

    def Compute(self, *args):
        return PyNormaliz_cpp.NmzCompute(self.cone, args)

    def setVerbose(self, verbose=True):
        return NmzSetVerbose(self.cone, verbose)

    # This one is not like the others!
    def IntegerHull(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IntegerHull")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        new_inner_cone = PyNormaliz_cpp.NmzResult(self.cone, "IntegerHull")
        return_cone = Cone.__new__(Cone)
        return_cone.cone = new_inner_cone
        return return_cone

    def ProjectCone(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ProjectCone")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        new_inner_cone = PyNormaliz_cpp.NmzResult(self.cone, "ProjectCone")
        return_cone = Cone.__new__(Cone)
        return_cone.cone = new_inner_cone
        return return_cone

    def EuclideanVolume(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Volume")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzGetEuclideanVolume(self.cone)

    def HilbertSeries(self, **kwargs):
        try:
            as_hsop = kwargs["HSOP"]
        except KeyError:
            as_hsop = 28
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("HilbertSeries")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        if as_hsop == 28:
            return PyNormaliz_cpp.NmzHilbertSeries(self.cone)
        if type(as_hsop) == bool:
            return PyNormaliz_cpp.NmzHilbertSeries(self.cone, as_hsop)
        raise TypeError("If HSOP is given, it must be True or False")

    def EhrhartSeries(self, **kwargs):
        try:
            as_hsop = kwargs["HSOP"]
        except KeyError:
            as_hsop = 28
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("EhrhartSeries")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        if as_hsop == 28:
            return PyNormaliz_cpp.NmzHilbertSeries(self.cone)
        if type(as_hsop) == bool:
            return PyNormaliz_cpp.NmzHilbertSeries(self.cone, as_hsop)
        raise TypeError("If HSOP is given, it must be True or False")

    def Polynomial(self, **kwargs):
        return PyNormaliz_cpp.NmzGetPolynomial(self.cone)

    def NrCoeffQuasiPol(self, bound=-1):
        return PyNormaliz_cpp.NmzSetNrCoeffQuasiPol(self.cone, bound)

    def SymmetrizedCone(self, **kwargs):
        new_inner_cone = PyNormaliz_cpp.NmzSymmetrizedCone(self.cone)
        if new_inner_cone == None:
            return None
        return_cone = Cone.__new__(Cone)
        return_cone.cone = new_inner_cone
        return return_cone

    def HilbertSeriesExpansion(self,degree):
        return NmzGetHilbertSeriesExpansion(self.cone,degree)

    def WeightedEhrhartSeriesExpansion(self,degree):
        return NmzGetWeightedEhrhartSeriesExpansion(self.cone,degree)

    def PrettyPolynomialTuple(self, numCoefficients, denCoefficients):
        """
        Strings for numerator and denominator of the a hilbert series.

        Parameters
        ----------
        numCoefficients : list
            The coefficients for the numerator.
        denCofficients : list
            The coefficients for the denominator where the value represents the
            exponent of 't' and the frequency indicates the outer coefficient.

        Returns
        -------
        PrettyPolynomialTuple: tuple of strings

        Examples
        --------

        >>> numCoefficients = [3, 7, 4, -4, -6, 5]
        >>> denCoefficients = [1, 1, 2, 2, 2, 4]
        >>> PrettyPolynomialTuple(numCoefficients,denCoefficients)

        ('(3 + 7t + 4t² - 4t³ - 6t⁴ + 5t⁵)', '(1 - t)² (1 - t²)³ (1 - t⁴)')

        """
        def to_sup(s):
            sups = {u'0': u'\u2070',
                    u'1': u'\xb9',
                    u'2': u'\xb2',
                    u'3': u'\xb3',
                    u'4': u'\u2074',
                    u'5': u'\u2075',
                    u'6': u'\u2076',
                    u'7': u'\u2077',
                    u'8': u'\u2078',
                    u'9': u'\u2079'}
            if s is 1:
                return ''
            # lose the list comprehension
            return ''.join(sups.get(str(char), str(char)) for char in str(s))

        def getNumerator(coefficients):

            numerator = ''

            def isPositive(x):
                return x > 0

            firstNonZero = next(
                (i for i, x in enumerate(coefficients) if x != 0), 0)
            for exp, coefficient in enumerate(coefficients):
                if coefficient is 0:
                    continue
                # Exponent is 0 so keep only the coefficient
                if exp is 0:
                    numerator += '({}{!s}'.format('-' if not isPositive(coefficient)
                                                  else '', abs(coefficient))
                # Only include sign if `coefficient` is negative
                elif i is firstNonZero:
                    numerator += '{}{!s}t{}'.format('-' if not isPositive(
                        coefficient) else '', abs(coefficient), to_sup(exp))
                else:
                    numerator += ' {}{!s}t{}'.format('+ ' if isPositive(
                        coefficient) else '- ', abs(coefficient), to_sup(exp))
            numerator += ')'
            return numerator

        def getDenominator(coefficients):
            exponents = [(inner, coefficients.count(inner))
                         for inner in set(coefficients)]
            denominator = ' '.join('(1 - t{}){}'. format(to_sup(x[0]) if x[
                                   0] is not 1 else '', to_sup(x[1]) if x[1] is not 1 else '') for x in exponents)
            return denominator

        num = getNumerator(numCoefficients)
        den = getDenominator(denCoefficients)
        prettyPolynomial = (num, den)
        return prettyPolynomial

    def PrintPrettyHilbertSeries(self, numCoefficients, denCoefficients):
        """
        Make a pretty hilbert series string

        Parameters
        ----------
        numCoefficients : list of ints
            The coefficients for the numerator.
        denCofficients : list of ints
            The coefficients for the denominator where the value represents
            the exponent of 't' and the frequency indicates the outer
            coefficient.

        Returns
        -------
        PrintPrettyHilbertSeries : string

        Examples
        --------

        >>> numCoefficients = [3, 7, 4, -4, -6, 5]
        >>> deCoefficients = [1, 1, 2, 2, 2, 4]
        >>> PrintPrettyHilbertSeries(numCoefficients,deCoefficients)

        (3 + 7t + 4t² - 4t³ - 6t⁴ + 5t⁵)
        --------------------------------
           (1 - t)² (1 - t²)³ (1 - t⁴)

        """
        num, den = self.PrettyPolynomialTuple(numCoefficients, denCoefficients)
        prettyPolynomial = '{:^}\n{:-^{width}}\n{:^{width}}'.format(
            num, '', den, width=max(len(den),len(num)))
        return prettyPolynomial

    def PrintHilbertSeries(self):
        hilbert_series=self.HilbertSeries()
        shift=hilbert_series[2]
        shift=[ 0 for x in range(1,shift) ]
        numerator=shift+hilbert_series[0]
        denominator=hilbert_series[1]
        print(self.PrintPrettyHilbertSeries(numerator,denominator))
        return None

    # Auto generated stuff

    def Generators(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Generators")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Generators")

    def ExtremeRays(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ExtremeRays")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ExtremeRays")

    def VerticesOfPolyhedron(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("VerticesOfPolyhedron")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "VerticesOfPolyhedron")

    def SupportHyperplanes(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("SupportHyperplanes")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "SupportHyperplanes")

    def HilbertBasis(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("HilbertBasis")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "HilbertBasis")

    def ModuleGenerators(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ModuleGenerators")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ModuleGenerators")

    def Deg1Elements(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Deg1Elements")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Deg1Elements")

    def ModuleGeneratorsOverOriginalMonoid(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ModuleGeneratorsOverOriginalMonoid")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ModuleGeneratorsOverOriginalMonoid")

    def Sublattice(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Sublattice")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Sublattice")

    def ExcludedFaces(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ExcludedFaces")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ExcludedFaces")

    def OriginalMonoidGenerators(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("OriginalMonoidGenerators")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "OriginalMonoidGenerators")

    def MaximalSubspace(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("MaximalSubspace")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "MaximalSubspace")

    def Equations(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Equations")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Equations")

    def Congruences(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Congruences")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Congruences")

    def Grading(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Grading")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Grading")

    def Dehomogenization(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Dehomogenization")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Dehomogenization")

    def WitnessNotIntegrallyClosed(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("WitnessNotIntegrallyClosed")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "WitnessNotIntegrallyClosed")

    def TriangulationSize(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("TriangulationSize")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "TriangulationSize")

    def TriangulationDetSum(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("TriangulationDetSum")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "TriangulationDetSum")

    def ReesPrimaryMultiplicity(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ReesPrimaryMultiplicity")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ReesPrimaryMultiplicity")

    def GradingDenom(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("GradingDenom")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "GradingDenom")

    def UnitGroupIndex(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("UnitGroupIndex")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "UnitGroupIndex")

    def InternalIndex(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("InternalIndex")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "InternalIndex")

    def ExternalIndex(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ExternalIndex")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ExternalIndex")

    def Multiplicity(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Multiplicity")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Multiplicity")

    def RecessionRank(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("RecessionRank")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "RecessionRank")

    def AffineDim(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("AffineDim")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "AffineDim")

    def ModuleRank(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ModuleRank")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ModuleRank")

    def Rank(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Rank")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Rank")

    def EmbeddingDim(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("EmbeddingDim")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "EmbeddingDim")

    def IsPointed(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsPointed")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsPointed")

    def IsDeg1ExtremeRays(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsDeg1ExtremeRays")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsDeg1ExtremeRays")

    def IsDeg1HilbertBasis(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsDeg1HilbertBasis")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsDeg1HilbertBasis")

    def IsIntegrallyClosed(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsIntegrallyClosed")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsIntegrallyClosed")

    def IsReesPrimary(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsReesPrimary")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsReesPrimary")

    def IsInhomogeneous(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsInhomogeneous")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsInhomogeneous")

    def Triangulation(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Triangulation")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Triangulation")

    def InclusionExclusionData(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("InclusionExclusionData")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "InclusionExclusionData")

    def StanleyDec(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("StanleyDec")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "StanleyDec")

    def ClassGroup(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ClassGroup")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ClassGroup")

    def ConeDecomposition(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("ConeDecomposition")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "ConeDecomposition")

    def HilbertQuasiPolynomial(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("HilbertQuasiPolynomial")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "HilbertQuasiPolynomial")

    def EhrhartQuasiPolynomial(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("HilbertQuasiPolynomial")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "HilbertQuasiPolynomial")

    def IsTriangulationNested(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsTriangulationNested")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsTriangulationNested")

    def IsTriangulationPartial(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsTriangulationPartial")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsTriangulationPartial")

    def WeightedEhrhartQuasiPolynomial(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("WeightedEhrhartQuasiPolynomial")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "WeightedEhrhartQuasiPolynomial")

    def WeightedEhrhartSeries(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("WeightedEhrhartSeries")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "WeightedEhrhartSeries")

    def Integral(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Integral")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Integral")

    def VirtualMultiplicity(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("VirtualMultiplicity")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "VirtualMultiplicity")

    def IsGorenstein(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("IsGorenstein")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "IsGorenstein")

    def GeneratorOfInterior(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("GeneratorOfInterior")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "GeneratorOfInterior")

    def VerticesFloat(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("VerticesFloat")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "VerticesFloat")

    def Volume(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("Volume")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "Volume")

    def SuppHypsFloat(self, **kwargs):
        input_list = self.__process_keyword_args(kwargs)
        input_list.append("SuppHypsFloat")
        PyNormaliz_cpp.NmzCompute(self.cone, input_list)
        return PyNormaliz_cpp.NmzResult(self.cone, "SuppHypsFloat")
