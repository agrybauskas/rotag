#ifndef _CIF_H_
#define _CIF_H_

extern "C" {
  #include <datablock.h>
  #include <cifmessage.h>
}

typedef struct CIF CIF;

struct CIF {
  int nerrors;
  int yyretval;
  int major_version;
  int minor_version;
  DATABLOCK *datablock_list;
  DATABLOCK *last_datablock; /* points to the end of the
                                datablock_list; SHOULD not be freed
                                when the CIF structure is deleted.*/

  DATABLOCK *current_datablock; /* points to the data block (which
                                   can represent a data block or a
                                   save frame) that is currently
                                   parsed; SHOULD not be freed when
                                   the CIF structure is deleted.*/

  CIFMESSAGE *messages; /* A linked list with error and warning
                           message data. */
};

#endif
