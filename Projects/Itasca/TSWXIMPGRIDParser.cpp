#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "pzerror.h"

#include "TSWXIMPGRIDParser.h"

using namespace std;

TSWXGridData * swx::ReadGridDataFile(const std::string &filename)
{
    time_t start = time(NULL);
    TSWXGridData * grid_data = new TSWXGridData();

    std::ifstream input_file(filename.c_str());
    std::vector< unsigned int > indices;
    indices.reserve(8);
    if (input_file.is_open())
    {
        string group_name;
        bool read_group = false;
        std::vector< unsigned int > elements;

        while(!input_file.eof())
        {
            string line;
            getline(input_file,line);
            stringstream stream(line);

            if (!read_group)
            {
                string flag;
                stream >> flag;

                if (flag == "G")
                {
                    unsigned int id;
                    double x, y, z;

                    stream >> id;
                    stream >> x;
                    stream >> y;
                    stream >> z;
                    grid_data->AddGridPoint(id, x, y, z);
                }
                else if (flag == "Z")
                {
                    int item;
                    string type;
                    unsigned int id;

                    stream >> type;
                    stream >> id;

                    indices.resize(0);
                    stream >> item;
                    indices.push_back(item);
                    stream >> item;
                    indices.push_back(item);
                    stream >> item;
                    indices.push_back(item);
                    stream >> item;
                    indices.push_back(item);
                    stream >> item;
                    indices.push_back(item);
                    stream >> item;
                    indices.push_back(item);

                    if (type == "B8")
                    {
                      stream >> item;
                      indices.push_back(item);
                      stream >> item;
                      indices.push_back(item);
                      grid_data->AddB8Element(id, indices);
                    }
                    else if (type == "W6")
                    {
                      indices.push_back(item);
                      grid_data->AddW6Element(id, indices);
                    }

                }
                else if (flag == "ZGROUP")
                {
                    stream >> group_name;
                    read_group = true;
                }
            }
            else
            {
                string item;
                stream >> item;
                
                if (item == "ZGROUP")
                {
                    grid_data->AddZoneGroup(group_name, elements);
                    stream >> group_name;
                    elements.clear();
                    read_group = true;
                }
                else
                {
                  if(!item.empty())
                  {
                    int value = atoi(item.c_str());
                    elements.push_back(value);

                    while(stream >> value)
                    {
                      elements.push_back(value);
                    }
                  }
                }
            }

        }   
        if (read_group)
        {
            grid_data->AddZoneGroup(group_name, elements);
            group_name = "";
            elements.clear();
            read_group = false;
        }
    }
    input_file.close();

   	time_t end = time(NULL);
    const double elapsed = difftime(end, start);

    cout << "Elapsed: " << elapsed << endl;
    cout << "Grid Point size: " << grid_data->GetSizeGridPoint() << endl;
    cout << "B8 Element size: " << grid_data->GetSizeB8Element() << endl;
    cout << "W6 Element size: " << grid_data->GetSizeW6Element() << endl;
    cout << "Zone Group size: " << grid_data->GetSizeZoneGroup() << endl;

    return grid_data;
}




void swx::WriteGridDataFile(const TSWXGridData * gridData, std::string &filename)
{
    std::ofstream output_file(filename.c_str());

    if(!gridData)
    {
      DebugStop();
    }

    output_file << "* FLAC3D grid produced by www.simworx.com.br\n";
    output_file << "* 2014\n";

    output_file << "* GRIDPOINTS\n";
    for(int p = 0; p < gridData->GetSizeGridPoint(); p++)
    {
      output_file << "G\t" <<
                     gridData->GridPoints()[p].fId   << "\t" <<
                     gridData->GridPoints()[p].fX[0]<< "\t"  <<
                     gridData->GridPoints()[p].fX[1]<< "\t"  <<
                     gridData->GridPoints()[p].fX[2]<< "\n";
    }

    output_file << "* ZONES\n";
    for(int z = 0; z < gridData->GetSizeB8Element(); z++)
    {
      output_file << "Z\tB8\t" <<
                     gridData->B8Elements()[z].fId   << "\t" <<
                     gridData->B8Elements()[z].fIndices[0]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[1]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[2]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[3]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[4]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[5]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[6]<< "\t"  <<
                     gridData->B8Elements()[z].fIndices[7]<< "\n";
    }
    for(int z = 0; z < gridData->GetSizeW6Element(); z++)
    {
      output_file << "Z\tW6\t" <<
                     gridData->W6Elements()[z].fId   << "\t" <<
                     gridData->W6Elements()[z].fIndices[0]<< "\t"  <<
                     gridData->W6Elements()[z].fIndices[1]<< "\t"  <<
                     gridData->W6Elements()[z].fIndices[2]<< "\t"  <<
                     gridData->W6Elements()[z].fIndices[3]<< "\t"  <<
                     gridData->W6Elements()[z].fIndices[4]<< "\t"  <<
                     gridData->W6Elements()[z].fIndices[5]<< "\n";
    }

    if(gridData->GetSizeZoneGroup() > 0)
    {
      output_file << "* GROUPS\n";
      for(int zg = 0; zg < gridData->GetSizeZoneGroup(); zg++)
      {
          output_file << "ZGROUP " << gridData->ZoneGroups()[zg].fGroupName << "\n";

          for(int elId = 0; elId < gridData->ZoneGroups()[zg].fEls.size(); elId++)
          {
              output_file << gridData->ZoneGroups()[zg].fEls[elId];
              if((elId+1) < 15 || (elId+1)%15 != 0) output_file << "\t";
              else output_file << "\n";
          }
      }
    }
    output_file.close();
}
