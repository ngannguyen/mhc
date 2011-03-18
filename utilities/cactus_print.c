#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
//#include "psl.h"

/*
 * The script outputs a maf file containing all the block in a flower and its descendants.
 */

void getFlowers(Flower *flower, FILE *fileHandle);
void getBlockInstances(Block *block, FILE *fileHandle);

char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if(strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) *(1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    }
    else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

void getMAFBlock(Block *block, FILE *fileHandle) {
    /*
     * Outputs a MAF representation of the block to the given file handle.
     */
    fprintf(fileHandle, "a score=%i\n", block_getLength(block) *block_getInstanceNumber(block));
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    while((segment = block_getNext(instanceIterator)) != NULL) {
        Sequence *sequence = segment_getSequence(segment);
        if(sequence != NULL) {
            char *sequenceHeader = formatSequenceHeader(sequence);
            int32_t start;
            if(segment_getStrand(segment)) {
                start = segment_getStart(segment) - sequence_getStart(sequence);
            }
            else { //start with respect to the start of the reverse complement sequence
                start = (sequence_getStart(sequence) + sequence_getLength(sequence) - 1) - segment_getStart(segment);
            }
            int32_t length = segment_getLength(segment);
            char *strand = segment_getStrand(segment) ? "+" : "-";
            int32_t sequenceLength = sequence_getLength(sequence);
            char *instanceString = segment_getString(segment);
            fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceHeader, start, length, strand, sequenceLength, instanceString);
            free(instanceString);
            free(sequenceHeader);
        }
    }
    block_destructInstanceIterator(instanceIterator);
}

void getGroupEnds(Group *group, FILE *fileHandle){
    Group_EndIterator *endIterator = group_getEndIterator(group);
    End * end;
    while((end = group_getNextEnd(endIterator))!= NULL){
        fprintf(fileHandle, "%s, ", cactusMisc_nameToString(end_getName(end)));
    }
    fprintf(fileHandle, "\n");
    group_destructEndIterator(endIterator);
    return;
}
void getChainLinks(Chain *chain, FILE *fileHandle){
    int32_t i;
    for (i=0; i<chain_getLength(chain); i++){
        Link * link = chain_getLink(chain,i);
        Group *linkgroup = link_getGroup(link);
        fprintf(fileHandle, "\t\tLink %d <%s - %s>: ", i, cactusMisc_nameToString(end_getName(link_get5End(link))), cactusMisc_nameToString(end_getName(link_get3End(link))));
        fprintf(fileHandle, "group %s\n", cactusMisc_nameToString(group_getName(linkgroup)));
    }
}
void getFlowerChains(Flower *flower, FILE *fileHandle){
    //CHAINS
    fprintf(fileHandle, "CHAINIterator:\n");
    Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
    Chain *chain;
    while((chain= flower_getNextChain(chainIterator)) != NULL){
        fprintf(fileHandle, "\tChain %s\n", cactusMisc_nameToString(chain_getName(chain)));
        int32_t j,k;
        //ChainLinks (groups)
        getChainLinks(chain, fileHandle);
        //ChainBlocks
        fprintf(fileHandle, "\t\tChainBlocks:\n");
            Block **blocks = chain_getBlockChain(chain, &j); 
        for(k=0; k<j; k++) {
           fprintf(fileHandle,"\t\t\tcBlock %s ",cactusMisc_nameToString(block_getName(blocks[k])));
           fprintf(fileHandle, "<%s - %s>\n",cactusMisc_nameToString(end_getName(block_get5End(blocks[k]))), cactusMisc_nameToString(end_getName(block_get3End(blocks[k]))));
           getBlockInstances(blocks[k],fileHandle);
           fprintf(fileHandle, "\n");
        }
        fprintf(fileHandle, "\n");
        free(blocks);
    }
    flower_destructChainIterator(chainIterator);
}
void getGroupLink(Group *group, FILE *fileHandle){
    Link * link = group_getLink(group);
    if (link != NULL){
        fprintf(fileHandle,"<%s - %s>",cactusMisc_nameToString(end_getName(link_get5End(link))), cactusMisc_nameToString(end_getName(link_get3End(link))));
    }
}
void getFlowerGroups(Flower *flower, FILE *fileHandle){
    //Call child flowers recursively.
    //GROUPS
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
        fprintf(fileHandle, "GROUPIterator:\n");
    while((group = flower_getNextGroup(groupIterator)) != NULL) {
                fprintf(fileHandle, "\tGroup %s; parent %s ", cactusMisc_nameToString(group_getName(group)), cactusMisc_nameToString(flower_getName(flower)));
        getGroupLink(group, fileHandle);
        fprintf(fileHandle, "\t");
                Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower != NULL) {
            fprintf(fileHandle, "EndIterator:\t");
            getGroupEnds(group, fileHandle);
                fprintf(fileHandle, "\t\thas nested flower. Enter nested flower...\n");
            getFlowers(group_getNestedFlower(group), fileHandle); //recursive call.
            fprintf(fileHandle, "\t\tEnd nested flower of group %s\n",cactusMisc_nameToString(group_getName(group)));
        }else{
                fprintf(fileHandle, "terminal, EndIterator:\t");
            getGroupEnds(group, fileHandle);
        }
    }
    flower_destructGroupIterator(groupIterator);
}

void getBlockInstances(Block *block, FILE *fileHandle){
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment * segment;
    //fprintf(fileHandle, "\t\t");
    Chain *chain = block_getChain(block);
    if(chain == NULL){
        fprintf(fileHandle, "Chain: NA\n");
    }else{
        fprintf(fileHandle, "Chain: %s\n", cactusMisc_nameToString(chain_getName(chain)));
    }
    while((segment = block_getNext(instanceIterator)) !=NULL){
        Sequence * sequence = segment_getSequence(segment);
        if (sequence != NULL){
            char *sequenceHeader = formatSequenceHeader(sequence);
            int length = segment_getLength(segment);
            int start = segment_getStart(segment);
            start = segment_getStrand(segment) ? start -2 : start - length;
            char strand = segment_getStrand(segment) ? '+' : '-'; 
            fprintf(fileHandle, "%s\t%d\t%d\t%c\t", sequenceHeader, start, length, strand);
            /*if (segment_get5Cap(segment) != NULL){
               fprintf(fileHandle, "5adj: %s-%d\t", cactusMisc_nameToString(cap_getName(cap_getAdjacency(segment_get5Cap(segment)))),
                                cap_getCoordinate(cap_getAdjacency(segment_get5Cap(segment))));
               fprintf(fileHandle, "5end: %s-%d\t", cactusMisc_nameToString(cap_getName(segment_get5Cap(segment))),
                                cap_getCoordinate(segment_get5Cap(segment)));
            }
            if (segment_get3Cap(segment) != NULL){
               fprintf(fileHandle, "3end: %s-%d\t", cactusMisc_nameToString(cap_getName(segment_get3Cap(segment))),
                                cap_getCoordinate(segment_get3Cap(segment)));
               fprintf(fileHandle, "3adj: %s-%d", cactusMisc_nameToString(cap_getName(cap_getAdjacency(segment_get3Cap(segment)))),
                                cap_getCoordinate(cap_getAdjacency(segment_get3Cap(segment))));
            }*/
            fprintf(fileHandle, "\n");
        } 
    }
    //fprintf(fileHandle, "\n");
    block_destructInstanceIterator(instanceIterator);
}

void getFlowerBlocks(Flower *flower, FILE *fileHandle){
        //BLOCKS
        fprintf(fileHandle, "BLOCKIterator:\n");
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while((block = flower_getNextBlock(blockIterator)) != NULL) {
                if(block_getChain(block) != NULL){
           fprintf(fileHandle, "\tc");
                }
        fprintf(fileHandle, "\tBlock %s:", cactusMisc_nameToString(block_getName(block)));
        fprintf(fileHandle, "<%s - %s> \n",cactusMisc_nameToString(end_getName(block_get5End(block))), cactusMisc_nameToString(end_getName(block_get3End(block))));
        //getBlockInstances(block, fileHandle);
        //getMAFBlock(block, fileHandle);
    }
    fprintf(fileHandle, "\n");
    flower_destructBlockIterator(blockIterator);
}

void getFlowerSequences(Flower *flower, FILE *fileHandle){
    //SEQUENCES
    fprintf(fileHandle, "SEQUENCEIterator:\n");
    Flower_SequenceIterator *sequenceIterator = flower_getSequenceIterator(flower);
    Sequence * sequence;
    while((sequence = flower_getNextSequence(sequenceIterator)) != NULL){
        //fprintf(fileHandle, "\tSequence %s, length %d\n", cactusMisc_nameToString(sequence_getName(sequence)), sequence_getLength(sequence));
        fprintf(fileHandle, "\tSequence %s, length %d\n", sequence_getHeader(sequence), sequence_getLength(sequence));
    }
    flower_destructSequenceIterator(sequenceIterator);
}

void getFlowerCaps(Flower *flower, FILE *fileHandle){
    fprintf(fileHandle, "CAPIterator: Name, strand, side, coordinate\n");
    Flower_CapIterator *instanceIterator = flower_getCapIterator(flower);
    Cap * instance;
    while((instance = flower_getNextCap(instanceIterator)) != NULL){
                if(end_isStubEnd(cap_getEnd(instance))){
                    fprintf(fileHandle, "STUB");
        } 
        char side = cap_getSide(instance) ? '5' : '3';
        char strand = cap_getStrand(instance) ? '+' : '-';
        fprintf(fileHandle,"\t%s %c %c %d\t",cactusMisc_nameToString(cap_getName(instance)), strand, side, cap_getCoordinate(instance));
        if(cap_getAdjacency(instance) == NULL){ continue; }
        fprintf(fileHandle,"adj: %s\n",cactusMisc_nameToString(cap_getName(cap_getAdjacency(instance))));
        Segment *segment = cap_getSegment(instance);
                if (segment != NULL){
                    int start = segment_getStart(segment);
                    int end = start + segment_getLength(segment);
                    fprintf(fileHandle, "segmentStart-End: %d\t%d\n",start, end);
        }
        //PRINTING CHILDREN
                /*int numChild = cap_getChildNumber(instance);
        int i = 0;
        for(i=0; i< numChild; i++){
                    Cap *child = cap_getChild(instance, i);
            if(child != NULL){
                char cside = cap_getSide(child) ? '5' : '3';
                char cstrand = cap_getStrand(child) ? '+' : '-';
                fprintf(fileHandle,"CHILD %d: %s %c %c %d\n", i, cactusMisc_nameToString(cap_getName(child)), cstrand, cside, cap_getCoordinate(child));
                  }
        }*/
    }
    fprintf(fileHandle, "\n");
    flower_destructCapIterator(instanceIterator);
}

void getFlowerEnds(Flower *flower, FILE *fileHandle){
    fprintf(fileHandle, "ENDIterator:\n");
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while((end = flower_getNextEnd(endIterator)) != NULL){
        fprintf(fileHandle,"%s",cactusMisc_nameToString(end_getName(end)));
        //if(end_isCap(end)){
        //    fprintf(fileHandle,"-CAP");
        //}
        if(end_isStubEnd(end)){
            fprintf(fileHandle,"-stub");
        }
        if(end_isBlockEnd(end)){
            fprintf(fileHandle,"-blockend");
        }
        fprintf(fileHandle,"\t");
    }
    fprintf(fileHandle, "\n");
    flower_destructEndIterator(endIterator);
}

End *getPseudoAdjacentEnd(End *end){
    assert(end != NULL);
    PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
    assert(pseudoAdjacency != NULL);
    assert(pseudoAdjacency_get3End(pseudoAdjacency) == end || pseudoAdjacency_get5End(pseudoAdjacency) == end);
    if (pseudoAdjacency_get3End(pseudoAdjacency) == end) {
        end = pseudoAdjacency_get5End(pseudoAdjacency);
    } else {
        end = pseudoAdjacency_get3End(pseudoAdjacency);
    }
    return end;
}

void getFlowerReference(Flower *flower, FILE *fileHandle){
    fprintf(fileHandle, "REFERENCE\n");
    Reference *reference = flower_getReference(flower);
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {
        fprintf(fileHandle, "\npseudoChrom %s\n", cactusMisc_nameToString(pseudoChromosome_getName(pseudoChromosome)) );
        End *end5 = pseudoChromosome_get5End(pseudoChromosome);
        //end5 = getPseudoAdjacentEnd(end5);
        End_InstanceIterator *it5 = end_getInstanceIterator(end5);
        Cap *cap5;
        fprintf(fileHandle, "\tStart:\n");
        while((cap5 = end_getNext(it5)) != NULL){
            Sequence *sequence5 = cap_getSequence(cap5);
            if(sequence5 != NULL){
                fprintf(fileHandle, "\t%s: %d\n", sequence_getHeader(sequence5), cap_getCoordinate(cap5));
            }
        }
        end_destructInstanceIterator(it5);

        End *end3 = pseudoChromosome_get3End(pseudoChromosome);
        //end3 = getPseudoAdjacentEnd(end3);
        End_InstanceIterator *it3 = end_getInstanceIterator(end3);
        Cap *cap3;
        fprintf(fileHandle, "\tEnd:\n");
        while((cap3 = end_getNext(it3)) != NULL){
            Sequence *sequence3 = cap_getSequence(cap3);
            if(sequence3 != NULL){
                fprintf(fileHandle, "\t%s: %d\n", sequence_getHeader(sequence3), cap_getCoordinate(cap3));
            }
        }

    }
    reference_destructPseudoChromosomeIterator(it);
}

void getFlowers(Flower *flower, FILE *fileHandle) {
        fprintf(fileHandle, "FLOWER %s:\n", cactusMisc_nameToString(flower_getName(flower)));
    getFlowerChains(flower, fileHandle);
    //getFlowerSequences(flower, fileHandle);
    //getFlowerCaps(flower, fileHandle);
    //getFlowerEnds(flower, fileHandle);
    //getFlowerBlocks(flower, fileHandle);
    //getFlowerGroups(flower, fileHandle);
    //getFlowerReference(flower, fileHandle);
}

void usage() {
    fprintf(stderr, "cactus_mafGenerator, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the MAFs in.\n");
    fprintf(stderr, "-f --includeTrees : Include trees for each MAF block inside of a comment line.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    CactusDisk *cactusDisk;
    Flower *flower;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "cactusDisk", required_argument, 0, 'c' },
            { "flowerName", required_argument, 0, 'd' },
            { "outputFile", required_argument, 0, 'e' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:h", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                flowerName = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
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

    assert(cactusDiskDatabaseString != NULL);
    assert(flowerName != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);
    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output MAF file : %s\n", outputFile);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////
    flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(flowerName));
    st_logInfo("Parsed the top level flower of the cactus tree to check\n");

    ///////////////////////////////////////////////////////////////////////////
    // Recursive check the flowers.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
        fprintf(stderr, "getting flower\n");
    //getFlowers(flower, stderr);
    flower_check(flower);

    getFlowers(flower, fileHandle);
    fclose(fileHandle);
    st_logInfo("Got the mafs in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);

    return 0;
}
