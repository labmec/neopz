#ifndef CYCLE_TIMER_H

#include<string.h>

extern "C" {
  __inline__ uint64_t rdtsc(void) {
    uint32_t lo, hi;
    __asm__ __volatile__ (      // serialize
    "xorl %%eax,%%eax \n        cpuid"
    ::: "%rax", "%rbx", "%rcx", "%rdx");
    /* We cannot use "=A", since this would use %rax on x86_64 and 
       return only the lower 32bits of the TSC */
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return (uint64_t) hi << 32 | lo;
  }
}

class CycleTimer
{
 public:
  CycleTimer() {};
  
  void start() {
    cycles = rdtsc();
  }

  uint64_t stop() {
    uint64_t aux = rdtsc();
    cycles = cycles - aux;
    return cycles;
  }

  uint64_t getCycles() {return cycles;}

  uint64_t getUnits() {return cycles;}
    
  std::string getTime() {
    std::string str = uint64_to_string(cycles); 
    str += " cycles"; 
    return str;
  }

  private:

  uint64_t cycles;
    
  std::string uint64_to_string(uint64_t value ) {
    std::ostringstream os;
    os << value;
    return os.str();
  }

};

#endif
