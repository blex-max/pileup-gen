#include <string>

#include <htslib/sam.h>


namespace readops {

  bam1_t * from_ref (std::string ref, size_t rpos, size_t read_len);

  void ins (bam1_t * b, size_t pos, std::string ins);
  void del (bam1_t * b, size_t pos, size_t len);
  void sub (bam1_t * b, size_t pos, char base);
  // void clip (bam1_t * b);

  // void assign_qual (bam1_t * b, QualModel qual);

}
