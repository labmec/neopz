
#ifndef TSWXIMPGRIDParserH
#define TSWXIMPGRIDParserH

#include "TSWXIMPGRIDData.h"
#include "string"

namespace swx{

TSWXGridData * ReadGridDataFile(const std::string &filename);

void WriteGridDataFile(const TSWXGridData * gridData, std::string &filename);

int GetSizeFile(std::string filename);

};

#endif

