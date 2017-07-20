#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

source $AL_DIR_TOOLS/alea.config


##############################################################################
#############   Module 0: interactive setting
##############################################################################


#-----------------------------------------------------------------------------
# Modifies the setting file
#-----------------------------------------------------------------------------
function modifySetting {
    local PARAM_ANALYSIS_TYPE=$1
    local PARAM_ALIGNER_TYPE=$2
    local PARAM_METHOD_TYPE=$3
    
    cp $AL_DIR_TOOLS/alea.setting $AL_DIR_TOOLS/.alea.setting.old
    
    case $PARAM_ANALYSIS_TYPE in
        1 ) # in the case of ChIP-seq
            cat $AL_DIR_TOOLS/alea.setting | sed -e "
                s/AL_ANALYSIS_TYPE=.*/AL_ANALYSIS_TYPE=$PARAM_ANALYSIS_TYPE/
                s/AL_ALIGNER_TYPE_CHIP=.*/AL_ALIGNER_TYPE_CHIP=$PARAM_ALIGNER_TYPE/
                s/AL_METHOD_TYPE=.*/AL_METHOD_TYPE=$PARAM_METHOD_TYPE/
            " > $AL_DIR_TOOLS/.alea.setting.tmp
            ;;
        
        2 ) # in the case of RNA-seq
            cat $AL_DIR_TOOLS/alea.setting | sed -e "
                s/AL_ANALYSIS_TYPE=.*/AL_ANALYSIS_TYPE=$PARAM_ANALYSIS_TYPE/
                s/AL_ALIGNER_TYPE_RNA=.*/AL_ALIGNER_TYPE_RNA=$PARAM_ALIGNER_TYPE/
                s/AL_METHOD_TYPE=.*/AL_METHOD_TYPE=$PARAM_METHOD_TYPE/
            " > $AL_DIR_TOOLS/.alea.setting.tmp
            ;;
        
        3) # in the case of BS-seq
            cat $AL_DIR_TOOLS/alea.setting | sed -e "
                s/AL_ANALYSIS_TYPE=.*/AL_ANALYSIS_TYPE=$PARAM_ANALYSIS_TYPE/
                s/AL_ALIGNER_TYPE_BS=.*/AL_ALIGNER_TYPE_BS=$PARAM_ALIGNER_TYPE/
                s/AL_METHOD_TYPE=.*/AL_METHOD_TYPE=$PARAM_METHOD_TYPE/
            " > $AL_DIR_TOOLS/.alea.setting.tmp
            ;;
        
        * )
            echo "Invalid analysis type"
            exit 1
            ;;
    esac
    
    if [ $AL_DEBUG = 0 ]
    then
        mv $AL_DIR_TOOLS/.alea.setting.tmp $AL_DIR_TOOLS/alea.setting
    else
        cp $AL_DIR_TOOLS/.alea.setting.tmp $AL_DIR_TOOLS/alea.setting
    fi
    
}

#-----------------------------------------------------------------------------
# Modifies the config file
#-----------------------------------------------------------------------------
function modifyConfig {
    local PARAM_ANALYSIS_TYPE=$1
    local PARAM_ALIGNER_TYPE=$2
    local PARAM_METHOD_TYPE=$3
    
    cp $AL_DIR_TOOLS/alea.config $AL_DIR_TOOLS/.alea.config.old
    
    case $PARAM_ANALYSIS_TYPE in
        1 ) # in the case of ChIP-seq
            case $PARAM_ALIGNER_TYPE in
                1 ) # in the case of BWA
                    cat $AL_DIR_TOOLS/alea.config | sed -e "
                        s/AL_USE_BWA=[01]/AL_USE_BWA=1/
                        s/AL_USE_BOWTIE1=[01]/AL_USE_BOWTIE1=0/
                        s/AL_USE_BOWTIE2=[01]/AL_USE_BOWTIE2=0/
                        s/AL_USE_BISMARK=[01]/AL_USE_BISMARK=0/
                        s/AL_USE_STAR=[01]/AL_USE_STAR=0/
                        s/AL_USE_TOPHAT2=[01]/AL_USE_TOPHAT2=0/
                    " > $AL_DIR_TOOLS/.alea.config.tmp
                    ;;
                
                2 ) # in the case of Bowtie1
                    cat $AL_DIR_TOOLS/alea.config | sed -e "
                        s/AL_USE_BWA=[01]/AL_USE_BWA=0/
                        s/AL_USE_BOWTIE1=[01]/AL_USE_BOWTIE1=1/
                        s/AL_USE_BOWTIE2=[01]/AL_USE_BOWTIE2=0/
                        s/AL_USE_BISMARK=[01]/AL_USE_BISMARK=0/
                        s/AL_USE_STAR=[01]/AL_USE_STAR=0/
                        s/AL_USE_TOPHAT2=[01]/AL_USE_TOPHAT2=0/
                    " > $AL_DIR_TOOLS/.alea.config.tmp
                    ;;
                
                3 ) # in the case of Bowtie2
                    cat $AL_DIR_TOOLS/alea.config | sed -e "
                        s/AL_USE_BWA=[01]/AL_USE_BWA=0/
                        s/AL_USE_BOWTIE1=[01]/AL_USE_BOWTIE1=0/
                        s/AL_USE_BOWTIE2=[01]/AL_USE_BOWTIE2=1/
                        s/AL_USE_BISMARK=[01]/AL_USE_BISMARK=0/
                        s/AL_USE_STAR=[01]/AL_USE_STAR=0/
                        s/AL_USE_TOPHAT2=[01]/AL_USE_TOPHAT2=0/
                    " > $AL_DIR_TOOLS/.alea.config.tmp
                    ;;
                
                * )
                    echo "Invalid aligner type"
                    exit 1
                    ;;
            esac
            ;;
        
        2 ) # in the case of RNA-seq
            case $PARAM_ALIGNER_TYPE in
                
                1 ) # in the case of STAR
                    cat $AL_DIR_TOOLS/alea.config | sed -e "
                        s/AL_USE_BWA=[01]/AL_USE_BWA=0/
                        s/AL_USE_BOWTIE1=[01]/AL_USE_BOWTIE1=0/
                        s/AL_USE_BOWTIE2=[01]/AL_USE_BOWTIE2=0/
                        s/AL_USE_BISMARK=[01]/AL_USE_BISMARK=0/
                        s/AL_USE_STAR=[01]/AL_USE_STAR=1/
                        s/AL_USE_TOPHAT2=[01]/AL_USE_TOPHAT2=0/
                    " > $AL_DIR_TOOLS/.alea.config.tmp
                    ;;
                
                2 ) # in the case of TOPHAT2
                    cat $AL_DIR_TOOLS/alea.config | sed -e "
                        s/AL_USE_BWA=[01]/AL_USE_BWA=0/
                        s/AL_USE_BOWTIE1=[01]/AL_USE_BOWTIE1=0/
                        s/AL_USE_BOWTIE2=[01]/AL_USE_BOWTIE2=0/
                        s/AL_USE_BISMARK=[01]/AL_USE_BISMARK=0/
                        s/AL_USE_STAR=[01]/AL_USE_STAR=0/
                        s/AL_USE_TOPHAT2=[01]/AL_USE_TOPHAT2=1/
                    " > $AL_DIR_TOOLS/.alea.config.tmp
                    ;;
                
                * )
                    echo "Invalid aligner type"
                    exit 1
                    ;;
            esac
            ;;
        
        3 ) # in the case of BS-seq
            case $PARAM_ALIGNER_TYPE in
                1 ) # in the case of Bismark
                    cat $AL_DIR_TOOLS/alea.config | sed -e "
                        s/AL_USE_BWA=[01]/AL_USE_BWA=0/
                        s/AL_USE_BOWTIE1=[01]/AL_USE_BOWTIE1=0/
                        s/AL_USE_BOWTIE2=[01]/AL_USE_BOWTIE2=0/
                        s/AL_USE_BISMARK=[01]/AL_USE_BISMARK=1/
                        s/AL_USE_STAR=[01]/AL_USE_STAR=0/
                        s/AL_USE_TOPHAT2=[01]/AL_USE_TOPHAT2=0/
                    " > $AL_DIR_TOOLS/.alea.config.tmp
                    ;;
                
                * )
                    echo "Invalid aligner type"
                    exit 1
                    ;;
            esac
            ;;
        
        * )
            echo "Invalid analysis type"
            exit 1
            ;;
    esac
    
    if [ $AL_DEBUG = 0 ]
    then
        mv $AL_DIR_TOOLS/.alea.config.tmp $AL_DIR_TOOLS/alea.config
    else
        cp $AL_DIR_TOOLS/.alea.config.tmp $AL_DIR_TOOLS/alea.config
        mv $AL_DIR_TOOLS/.alea.config.tmp $AL_DIR_TOOLS/.alea.config.tmp1
    fi
    
    case $PARAM_METHOD_TYPE in
        1 ) # in the case of concatenated method
            cat $AL_DIR_TOOLS/alea.config | sed -e "
                s/AL_USE_CONCATENATED_GENOME=[01]/AL_USE_CONCATENATED_GENOME=1/
            " > $AL_DIR_TOOLS/.alea.config.tmp
            ;;
        
        2 ) # in the case of separate method
            cat $AL_DIR_TOOLS/alea.config | sed -e "
                s/AL_USE_CONCATENATED_GENOME=[01]/AL_USE_CONCATENATED_GENOME=0/
            " > $AL_DIR_TOOLS/.alea.config.tmp
            ;;
        
        * )
            echo "Invalid method type"
            exit 1
            ;;
    esac
    
    if [ $AL_DEBUG = 0 ]
    then
        #mv $AL_DIR_TOOLS/.alea.config.tmp $AL_DIR_TOOLS/alea.config
        mv $AL_DIR_TOOLS/.alea.config.tmp /alea-data/alea.config
    else
        cp $AL_DIR_TOOLS/.alea.config.tmp $AL_DIR_TOOLS/alea.config
        mv $AL_DIR_TOOLS/.alea.config.tmp $AL_DIR_TOOLS/.alea.config.tmp2
    fi
}

#------------------------------------------------------------------------------------
# setting
#------------------------------------------------------------------------------------

echo
echo -n "Which type of analysis do you want to do?"
echo "$AL_ANALYSIS_CHOICE"
echo -n "Press 1, 2 or 3 [current: $AL_ANALYSIS_TYPE] : "

read INPUT_ANALYSIS_TYPE
if [ -z "$INPUT_ANALYSIS_TYPE" ]
then
    INPUT_ANALYSIS_TYPE=$AL_ANALYSIS_TYPE
elif [[ ! "$INPUT_ANALYSIS_TYPE" =~ [1-3] ]]
then
    echo "Invalid input. The current status was held."
    INPUT_ANALYSIS_TYPE=$AL_ANALYSIS_TYPE
fi

case $INPUT_ANALYSIS_TYPE in
    1 ) # in the case of ChIP-seq
        echo 
        echo -n "Which type of alignment tool do you want to use?"
        echo "$AL_ALIGNER_CHOICE_CHIP"
        echo -n "Press 1, 2 or 3 [current: $AL_ALIGNER_TYPE_CHIP] : "
        
        read INPUT_ALIGNER_TYPE
        if [ -z "$INPUT_ALIGNER_TYPE" ]
        then
            INPUT_ALIGNER_TYPE=$AL_ALIGNER_TYPE_CHIP
        elif [[ ! "$INPUT_ALIGNER_TYPE" =~ [1-3] ]]
        then
            echo "Invalid input. The current status was held."
            INPUT_ALIGNER_TYPE=$AL_ALIGNER_TYPE_CHIP
        fi
        ;;
    
    2 ) # in the case of RNA-seq
        echo 
        echo -n "Which type of alignment tool do you want to use?"
        echo "$AL_ALIGNER_CHOICE_RNA"
        echo -n "Press 1, 2 or 3 [current: $AL_ALIGNER_TYPE_RNA] : "
        
        read INPUT_ALIGNER_TYPE
        if [ -z "$INPUT_ALIGNER_TYPE" ]
        then
            INPUT_ALIGNER_TYPE=$AL_ALIGNER_TYPE_RNA
        elif [[ ! "$INPUT_ALIGNER_TYPE" =~ [1-5] ]]
        then
            echo "Invalid input. The current status was held."
            INPUT_ALIGNER_TYPE=$AL_ALIGNER_TYPE_RNA
        fi
        ;;
    
    3 ) # in the case of BS-seq
        echo 
        echo -n "Which type of alignment tool do you want to use?"
        echo "$AL_ALIGNER_CHOICE_BS"
        echo -n "Press 1 [current: $AL_ALIGNER_TYPE_BS] : "
        
        read INPUT_ALIGNER_TYPE
        if [ -z "$INPUT_ALIGNER_TYPE" ]
        then
            INPUT_ALIGNER_TYPE=$AL_ALIGNER_TYPE_BS
        elif [[ ! "$INPUT_ALIGNER_TYPE" =~ [1] ]]
        then
            echo "Invalid input. The current status was held."
            INPUT_ALIGNER_TYPE=$AL_ALIGNER_TYPE_BS
        fi
        ;;
    
    * )
        echo "Invalid analysis type"
        exit 1
        ;;
esac

echo
echo -n "Which type of alignment method do you want to use?"
echo "$AL_METHOD_CHOICE"
echo -n "Press 1 or 2 [current: $AL_METHOD_TYPE] : "

read INPUT_METHOD_TYPE
if [ -z "$INPUT_METHOD_TYPE" ]
then
    INPUT_METHOD_TYPE=$AL_METHOD_TYPE
elif [[ ! "$INPUT_METHOD_TYPE" =~ [1-3] ]]
then
    echo "Invalid input. The current status was held."
    INPUT_METHOD_TYPE=$AL_METHOD_TYPE
fi

modifySetting "$INPUT_ANALYSIS_TYPE" "$INPUT_ALIGNER_TYPE" "$INPUT_METHOD_TYPE"
modifyConfig "$INPUT_ANALYSIS_TYPE" "$INPUT_ALIGNER_TYPE" "$INPUT_METHOD_TYPE"
