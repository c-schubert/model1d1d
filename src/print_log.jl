function print_log(logfile::String, printargs; write_log=true, print_to_stdout = true)

    if print_to_stdout
        println(printargs...)
    end

    if write_log
        open(logfile, "a") do lf
            # modify mile content
            for arg in printargs
                write_lf_entries(lf, arg)
            end
            write(lf, "\n")
        end
    end

    return nothing
end

function write_lf_entries(lf, arg)
    if typeof(arg) <: Array
        write(lf,"[")
        for (i,arg_i) in enumerate(arg)
            write_lf_entries(lf, arg_i) 
            if i != length(arg)
                write(lf,",")
            end
        end
        write(lf,"]")
    else
        if typeof(arg) <: Number
            write(lf,string(arg))
        elseif typeof(arg) == String
            write(lf, arg)
        elseif typeof(arg) == Char
            write(lf, string(arg))
        else
            println(typeof(arg))
            error("Unknow type for writing to log file")
        end
    end
    return nothing
end

function new_log(logfile::String)
    if isfile(logfile)
        println("Logfile \"",logfile ,"\"allready exist. Continue and overwrite? y/N");
        s = readline()
        if (s == "Y" || s ==  "y")
            rm(logfile)
        else
            error("Manual abort, due to existing logfile!")
        end
    end

    println("Set m1d1d_logfile varible to ", logfile)
    global m1d1d_logfile = logfile

    print_log(m1d1d_logfile, string(Dates.now()))
end

function test_global_logfilename()
    println(m1d1d_logfile)
end