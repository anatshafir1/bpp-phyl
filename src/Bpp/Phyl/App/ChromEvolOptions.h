#ifndef _CHROMEVOLOPTIONS_H_
#define _CHROMEVOLOPTIONS_H_

// From bpp-core:

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include "PhylogeneticsApplicationTools.h"



//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
namespace bpp{

class ChromEvolOptions
{

public:
    static void initAllParameters(BppApplication& ChromEvol);
    virtual ~ChromEvolOptions();
    
public:
    static TreeTemplate<Node>* tree_;
    static const VectorSiteContainer* vsc_;
    static const ChromosomeAlphabet* alpha_;
    static int maxChrNum_;
    static int minChrNum_;
    static double branchMul_;
    static std::vector <unsigned int> OptPointsNum_;
    static std::vector <unsigned int> OptIterNum_;
    static double constGain_;
    static double constLoss_;
    static double constDupl_;
    static double constDemiDupl_;
    static double tolerance_;
    static unsigned int maxIterations_;
    static bool maxParsimonyBound_;
    static unsigned int maxAlpha_;
    static unsigned int minAlpha_;

private:
    static void initDefaultParameters();
    static void initParametersFromFile(BppApplication& ChromEvol);
    static VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc);
    static VectorSiteContainer* getCharacterData(const std :: string &path, unsigned int* numberOfUniqueStates);
    static void setMaxChrNum(unsigned int maxNumberOfChr);
    static void setMinChrNum(unsigned int minNumberOfChr);
    static void rescale_tree(TreeTemplate<Node>* tree, unsigned int chrRange);
    static TreeTemplate<Node>* getTree(const std :: string &path, unsigned int numOfUniqueCharacterStates);

};

}


#endif  // _CHROMEVOLOPTIONS_H_