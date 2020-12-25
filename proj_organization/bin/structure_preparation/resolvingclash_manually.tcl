proc clashes {selection1 selection2} {
        set x {}
        foreach coord1 [$selection1 get {x y z}] name1 [$selection1 get name] rad1 [$selection1 get radius] {
                foreach coord2 [$selection2 get {x y z}] name2 [$selection2 get name] rad2 [$selection2 get radius] {
                        set atmdist [veclength [vecsub $coord1 $coord2]]
                        set atmrad [expr $rad1 + $rad2 - 0.5]
                        if { $atmrad - $atmdist > 0 } {
                                lappend x $atmrad - $atmdist $name1 $name2
                        }
                }
        }
        puts "test2"
        puts [llength $x]
        if {[llength $x] > 0} {
        puts [llength $x]
        }
}




set c [atomselect top "chain 'K' 'L' 'M' 'N' 'O' 'T' 'U'"]
set a [atomselect top "chain 'A' 'B' 'C' 'D' 'E' 'P' 'Q'"]
set b [atomselect top "chain 'F' 'G' 'H' 'J' 'I' 'R' 'S'"]

clashes $a $b
#830
$a moveby {1 0 0}
$a moveby {-1 0 0}
$b moveby {-1 0 0}
$b moveby {1 0 0}
$b moveby {1 0 0}
clashes $a $b

$b moveby {1 0 0}
clashes $a $b

#above moveby command can also be used to move strands and create nanocrystal with varied number of strands geometry on different places,
