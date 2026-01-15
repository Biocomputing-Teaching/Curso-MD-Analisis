#!/usr/bin/env vmd

set coursedir ""
if {[llength $argv] > 0} {
    set coursedir [file normalize [lindex $argv 0]]
} elseif {[info exists ::env(COURSE_DIR)]} {
    set coursedir [file normalize $::env(COURSE_DIR)]
} else {
    set coursedir [file normalize [file join $env(HOME) Concepcion26]]
}

set data_dir [file join $coursedir data complex]
set pdb_file [file join $data_dir protein2.pdb]
set traj_file [file join $coursedir results 03-simulaciones-clasicas complex output_traj.dcd]
set out_dir [file join $coursedir results 04-analisis-trayectorias complex vmd]
file mkdir -p $out_dir
set rmsd_file [file join $out_dir rmsd_vmd.csv]
set rmsf_file [file join $out_dir rmsf_vmd.csv]
set dist_file [file join $out_dir dist_vmd.csv]
set movie_dir [file join $out_dir movie]
file mkdir -p $movie_dir

mol new $pdb_file
mol addfile $traj_file waitfor all
set proteinCA [atomselect top "protein and name CA"]
set refCA [atomselect top "protein and name CA" frame 0]
set ligand_sel [atomselect top "resname SUB"]
set anchor_sel [atomselect top "protein and resid 5 and name CA"]
set our_ligand [atomselect top "resname SUB and not name H*"]

set nframes [molinfo top get numframes]
set dt_ps 0.002

set fout [open $rmsd_file w]
puts $fout "Time (ps),Backbone RMSD (Å)"
set dist_out [open $dist_file w]
puts $dist_out "Time (ps),Anchor-Ligand distance (Å)"
for {set i 0} {$i < $nframes} {incr i} {
    animate goto $i
    $proteinCA frame $i
    $refCA frame 0
    $anchor_sel frame $i
    $our_ligand frame $i
    $ligand_sel frame $i
    set rmsd [measure rmsd $proteinCA $refCA]
    set time [expr {$i * $dt_ps}]
    puts $fout "[format %.3f $time],[format %.4f $rmsd]"
    set prot_pos [$anchor_sel get {x y z}]
    set lig_pos [$our_ligand get {x y z}]
    set delta [vecsub $prot_pos $lig_pos]
    set dist [vecnorm $delta]
    set fname [file join $movie_dir [format "frame_%04d.png" $i]]
    render TachyonInternal $fname
}
close $fout
close $dist_out

set rmsf_data [measure rmsf $proteinCA]
set resids [$proteinCA get resid]
array set atommap {}
for {set idx 0} {$idx < [llength $resids]} {incr idx} {
    set resid [lindex $resids $idx]
    lappend atommap($resid) $idx
}

set fout_rmsf [open $rmsf_file w]
puts $fout_rmsf "Residue ID,RMSF (Å)"
foreach resid [lsort -integer -unique $resids] {
    set idxs $atommap($resid)
    set sum 0
    foreach atom_idx $idxs {
        set sum [expr {$sum + [lindex $rmsf_data $atom_idx]}]
    }
    set avg [expr {$sum / [llength $idxs]}]
    puts $fout_rmsf "$resid,[format %.4f $avg]"
}
close $fout_rmsf

puts "Saved RMSD to $rmsd_file"
puts "Saved RMSF to $rmsf_file"
puts "Saved distances to $dist_file"
puts "Rendered $nframes frames in $movie_dir"
