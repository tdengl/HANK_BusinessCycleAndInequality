#---------------------------------------------------------
# Pre-process aggregate-model and steady-state inputs and write 
# them into the respective functions.
#---------------------------------------------------------

# aggregate steady state
f = open("8_Rank/prepare_linearization_rank.jl")
lines = readlines(f)
insert_index = findall( x -> x == "    # aggregate steady state marker", lines)[1]

g = open("8_Rank/input_aggregate_steady_state_rank.jl")
s = read(g, String)

open("8_Rank/prepare_linearization_generated_rank.jl", "w") do h
    println(h, "# This file has been automatically generated by PreprocessInputs_rank.jl. Any user inputs might be overwritten!")
    println(h, "\n")
    
    for i = 1:insert_index-1
        println(h, lines[i])
    end
    println(h, "\n")
    write(h, s)
    println(h, "\n")
    for i = insert_index+1:length(lines)
        println(h, lines[i])
    end
end
close(f)
close(g)

# aggregate model
f = open("5_LinearizationFunctions/FSYS_agg.jl")
lines = readlines(f)
insert_index = findall( x -> x == "    # aggregate model marker", lines)[1]

g = open("8_Rank/input_aggregate_model_rank.jl")
s = read(g, String)

open("8_Rank/FSYS_agg_generated.jl", "w") do h
    println(h, "# This file has been automatically generated by PreprocessInputs.jl. Any user inputs might be overwritten!")
    println(h, "\n")

    for i = 1:insert_index-1
        println(h, lines[i])
    end
    println(h, "\n")
    write(h, s)
    println(h, "\n")
    for i = insert_index+1:length(lines)
        println(h, lines[i])
    end
end
close(f)
close(g)