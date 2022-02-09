
process IN_TEST4 {

    echo true
   
    """
    if [ ${params.emptyrun} ]
    then
        echo IN_TEST4 module OFF
    else
        echo IN_TEST4 module task ${params.dryrun}
    fi
    """
}
