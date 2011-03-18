#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "common.h"
#include "psl.h"

/*
 *This script converts genome coordinates of input PSLs (from a file) to cactus' coordinates.
 */

void usage();
//==============MAP LOCATION OF PSL BLOCKS TO THE REFERENCE SEQUENCES================
int psl_getStrand(char strand){
   if(strand == '-'){
      return 0;
   }else{
      return 1;
   }
}
/*
int32_t mapCoor(int32_t refsize, int32_t start, int32_t coor, int strand){
   //Map 'coor' to the coordinates of the reference sequence
    //refsize: size of reference sequence, start: location of current seq on the ref-seq
    //strand: 1 if current seq is on the forward strand of the ref-seq, 0 otherwise
    //
   if(strand == 0){
      return refsize -1 - (coor + start);
   }else{
      return coor + start;
   }
   //return strand? coor + start : refsize - (coor + start);   
}
int32_t mapStrand(bool refstrand, bool strand){
   return strand ? refstrand : !(refstrand);
}
*/

//Given that all of cactus' input sequences are on the forward strand
void mapPslCoor(struct psl *pslRecord, int32_t trefsize, int32_t tstart, int32_t tlen){
   int i;
   pslRecord->tSize = trefsize;
   pslRecord->tStart = pslRecord->tStart - tstart + 2;
   pslRecord->tEnd = pslRecord->tEnd - tstart + 2;
   for(i=0; i< pslRecord->blockCount; i++){
      *(pslRecord->tStarts + i) = *(pslRecord->tStarts +i) - tstart + 2;
   }
   return;
}

void mapPSLs(struct psl *psl, FILE *fileHandle, char *target) {
   /*
    *Map coordinates of all the psls in the list psl to the query and target reference genomes
    */
   int32_t trefsize = 0;
   int32_t tstart = 0;
   int32_t tlen = 0;
   char sep = '\t';
   char lastStep = '\n';
   if( sscanf(target, "%d.%d.%d", &trefsize, &tstart, &tlen) != 3){
      fprintf(stderr, "Wrong target format\n");
      usage();
      return;
   }
   while(psl != NULL){
      mapPslCoor(psl, trefsize, tstart, tlen);
      pslOutput(psl, fileHandle, sep, lastStep);
      psl = psl->next;
   }
   return;
}
void usage() {
   fprintf(stderr, "pslConvertCoor, version 0.0\n");
   fprintf(stderr, "Convert genome coordinates of input psls to cactus' coordinates.\n");
   fprintf(stderr, "Note: psls must come from sequences that are on the forward strand of the genome sequences.\n");
   fprintf(stderr, "-a --logLevel : Set the log level\n");
   fprintf(stderr, "-i --inputFile : The file that contains list of input PSLs.\n");
   fprintf(stderr, "-o --outputFile : The file to write the new PSLs in.\n");
   fprintf(stderr, "-t --target location on its reference genome (forward strand): chrSize.qStart.qLength. E.g: 171115067.31321610.3479\n");
   fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
   /*
    * Arguments/options
    */
   char outputFile[50];
   char inputFile[50];
   char target[100];

   ///////////////////////////////////////////////////////////////////////////
   // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
   ///////////////////////////////////////////////////////////////////////////

   while(1) {
      static struct option long_options[] = {
         { "target", required_argument, 0, 't' },
         { "outputFile", required_argument, 0, 'o' },
         { "inputFile", required_argument, 0, 'i' },
         { "help", no_argument, 0, 'h' },
         { 0, 0, 0, 0 }
      };

      int option_index = 0;

      int key = getopt_long(argc, argv, "i:o:t:h", long_options, &option_index);

      if(key == -1) {
         break;
      }

      switch(key) {
         case 'i':
            strcpy(inputFile, optarg);
            break;
         case 'o':
            strcpy(outputFile, optarg);
            break;
         case 't':
            strcpy(target, optarg);
            break;
         case 'h':
            usage();
            return 0;
         default:
            usage();
            return 1;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   // (0) Check the inputs.
   ///////////////////////////////////////////////////////////////////////////

   assert(outputFile != NULL);
   assert(target != NULL);

   FILE *fileHandle = fopen(outputFile, "w");
   pslWriteHead(fileHandle);
   struct psl *pslList = pslLoadAll(inputFile);
   mapPSLs(pslList, fileHandle, target);
   fclose(fileHandle);
   
   return 0;
}
