#ifndef TPZENUMAPPROXFAMILY_H
#define TPZENUMAPPROXFAMILY_H

/// Enum stating which flavor of HDiv spaces is being used
enum class HDivFamily {EHDivStandard,EHDivConstant,EHDivKernel};

/// Enum stating which flavor of H1 spaces is being used
enum class H1Family {EH1Standard};

/// Enum stating which flavor of HCurl spaces is being used
enum class HCurlFamily {EHCurlStandard,EHCurlNoGrads};

struct DefaultFamily {
    static const HDivFamily fHDivDefaultValue = HDivFamily::EHDivStandard;
    static const H1Family fH1DefaultValue = H1Family::EH1Standard;
    static const HCurlFamily fHCurlDefaultValue = HCurlFamily::EHCurlStandard;
};


#endif
