#include "MonoHiggsToGG/RooUtils/interface/RooPowLogPdf.h"
#include "MonoHiggsToGG/RooUtils/interface/RooSlicePdf.h"
#include "MonoHiggsToGG/RooUtils/interface/RooStarMomentMorph.h"

namespace  {
    struct dictionary {
	    RooPowLogPdf pl;
	    RooSlicePdf  sl;
	    RooStarMomentMorph  sm;
    };
}

