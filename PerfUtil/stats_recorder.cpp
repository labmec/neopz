#include "stats_recorder.h"
#include "pzerror.h"
#include <cstring> //memset
#ifdef USING_PAPI
#include <papi.h>
#endif

PAPIRunStat::PAPIRunStat() {
#ifndef USING_PAPI
  PZError << "For using PAPIRunStat use CMake option";
  PZError << " USING_PAPI" << std::endl;
  DebugStop();
#else
  clearStats();
#endif
}

void PAPIRunStat::start() {
#ifdef USING_PAPI
  PAPI_flops(&rtimeS, &ptimeS, &flpopsS, &mflopsS);
#endif
}

void PAPIRunStat::stop() {
#ifdef USING_PAPI
  float rtimeP, ptimeP, mflopsP;
  int64_t flpopsP;
  PAPI_flops(&rtimeP, &ptimeP, &flpopsP, &mflopsP);
  rtimeACC += (rtimeP - rtimeS);
  ptimeACC += (ptimeP - ptimeS);
  flpopsACC += (flpopsP - flpopsS);
#endif
}

void PAPIRunStat::print(std::ostream &os) const {
#ifdef USING_PAPI
  os << "HEADERS,PAPI_RTIME,PAPI_PTIME,PAPI_FLPOPS,PAPI_MFLOPS" << std::endl;
  os << "VALUES," << rtimeACC << "," << ptimeACC << "," << flpopsACC << ","
     << ((double)flpopsACC / (double)ptimeACC) / 1000000.0 << std::endl;
#endif
}
int PAPIRunStat::setCellValues(CSVStringTable &st, unsigned row) const {
#ifdef USING_PAPI
  if (st.nRows() <= row)
    return -1;
  int ret;
  if ((ret = st.setCell(row, "PAPI_RTIME", rtimeACC, true)))
    return ret;
  if ((ret = st.setCell(row, "PAPI_PTIME", ptimeACC, true)))
    return ret;
  if ((ret = st.setCell(row, "PAPI_FLPOPS", flpopsACC, true)))
    return ret;
  double av_mflops = ((double)flpopsACC / (double)ptimeACC) / 1000000.0;
  if ((ret = st.setCell(row, "PAPI_AV_MFLOPS", av_mflops, true)))
    return ret;
#endif
  return 0; // Return ok
}

void PAPIRunStat::clearStats() {
#ifdef USING_PAPI
  memset(&rtimeACC, 0, sizeof(rtimeACC));
  memset(&ptimeACC, 0, sizeof(ptimeACC));
  memset(&flpopsACC, 0, sizeof(flpopsACC));
#endif
}

ElapsedTimeRunStat::ElapsedTimeRunStat(){
  clearStats();
}

void ElapsedTimeRunStat::start() {
  gettimeofday(&lap, NULL);
}

void ElapsedTimeRunStat::stop() {
  timeval curr;
  gettimeofday(&curr, NULL);
  total.tv_sec += (curr.tv_sec - lap.tv_sec);
  total.tv_usec += (curr.tv_usec - lap.tv_usec);
  // elapsed = (t.tv_sec - lap.tv_sec) * 1000.0;      // sec to ms
  // elapsed += (t.tv_usec - lap.tv_usec) / 1000.0;   // us to ms
}

void ElapsedTimeRunStat::print(std::ostream &os) const {
  os << "HEADERS,ELAPSED_MS" << std::endl;
  os << "VALUES," << getElapsedMS() << std::endl;
}

int ElapsedTimeRunStat::setCellValues(CSVStringTable &st, unsigned row) const {
  if (st.nRows() <= row)
    return -1;
  return st.setCell(row, "ELAPSED", getElapsedMS(), true);
}

void ElapsedTimeRunStat::clearStats() {
  memset(&lap, 0, sizeof(lap));
  memset(&total, 0, sizeof(total));
}

double ElapsedTimeRunStat::getElapsedMS() const {
  double elapsed;
  elapsed = (total.tv_sec) * 1000.0;   // sec to ms
  elapsed += (total.tv_usec) / 1000.0; // us to ms
  return elapsed;
}

RunStatsRecorder::RunStatsRecorder() : n_laps(0) {
#ifdef USING_PAPI
  /* Add the PAPI counters statistics. */
  stat_items.push_back(new PAPIRunStat());
#endif
  /* Add the elapsed time statistic. */
  stat_items.push_back(new ElapsedTimeRunStat());
}

RunStatsRecorder::~RunStatsRecorder() {
  std::vector<RunStat *>::iterator it;
  for (it = stat_items.begin(); it != stat_items.end(); it++) {
    RunStat *i = *it;
    delete i;
  }
  /** @brief vector::clear Removes all elements from the vector (which are
   * destroyed), leaving the container with a size of 0. */
  stat_items.clear();
}

void RunStatsRecorder::start(){
  std::vector<RunStat *>::iterator it;
  for (it = stat_items.begin(); it != stat_items.end(); it++)
    (*it)->start();
}

void RunStatsRecorder::stop(){
  std::vector<RunStat *>::reverse_iterator rit;
  for (rit = stat_items.rbegin(); rit != stat_items.rend(); rit++)
    (*rit)->stop();
  n_laps++;
}

void RunStatsRecorder::lap(){
  stop();
  start();
}

void RunStatsRecorder::print(std::ostream& os) const{
  std::vector<RunStat*>::const_iterator it;
  for (it=stat_items.begin(); it!=stat_items.end(); it++)
    (*it)->print(os);
}

int RunStatsRecorder::append_to(CSVStringTable& st) const{
  int ret;
  unsigned new_row = st.addRows(1);

  if ((ret = update_row(st, new_row)))
    return ret;
  else
    return new_row;
}

int RunStatsRecorder::update_row(CSVStringTable &st, unsigned row) const {
  int ret;

  std::vector<RunStat *>::const_iterator it;
  for (it = stat_items.begin(); it != stat_items.end(); it++) {
    if ((ret = (*it)->setCellValues(st, row))) {
      return ret;
    }
  }

  return 0;
}

void RunStatsRecorder::clear() {
  std::vector<RunStat *>::iterator it;
  for (it = stat_items.begin(); it != stat_items.end(); it++)
    (*it)->clearStats();
  n_laps = 0;
}
