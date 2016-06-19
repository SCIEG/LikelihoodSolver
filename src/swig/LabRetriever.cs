//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.10
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------


public class LabRetriever {
  public static double CalculateLikelihood(bool calculateSuspect, uint numUnknowns, string alleleFrequencyTableDirectory, string locusName, string race, StringVector suspectAlleles, StringVectorVector assumedAlleles, StringVectorVector unattributedAlleles, double zeroAllelesInCommonProb, double oneAlleleInCommonProb, double bothAllelesInCommonProb, double dropoutRate, double dropinRate, double alpha, double fst) {
    double ret = LabRetrieverPINVOKE.CalculateLikelihood(calculateSuspect, numUnknowns, alleleFrequencyTableDirectory, locusName, race, StringVector.getCPtr(suspectAlleles), StringVectorVector.getCPtr(assumedAlleles), StringVectorVector.getCPtr(unattributedAlleles), zeroAllelesInCommonProb, oneAlleleInCommonProb, bothAllelesInCommonProb, dropoutRate, dropinRate, alpha, fst);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

}
