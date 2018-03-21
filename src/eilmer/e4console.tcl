#!/usr/bin/wish
#
# e4console.tcl
# GUI front-end for Eilmer4.
# P. A. Jacobs
#
# Versions:
# 04-Oct-2005 : First hack adapted from plot/sourse/mb_view.tcl
# 20-Jun-2006 : brought over to mbcns2
# 03-Apr-2008 : ported to Eilmer3
# 15-Aug-2008 : cleaned up to leave a minimal set of options.
# 05-Oct-2009 : Updated for the *current* Eilmer3 code.
# 19-Mar-2018 : Start refresh for Eilmer4.
#
# -----------------------------------------------------------------------

proc showAboutBox {} {
    tk_messageBox -type ok -title "About e4console.tcl" \
	-message "e4console V0.2\nP.J. 21-Mar-2018"
}; # end proc showAboutBox

proc showGeneralHelp {} {
    tk_messageBox -type ok -title "General Help" -message "None available yet."
    # Should change this procedure to insert text straight into the text widget.
}; # end proc showGeneralHelp

# -----------------------------------------------------------------------

# When looking at the OS name, we are expecting values such as
# "GNU/Linux" or "Cygwin"
if { [catch {exec uname -o} ::my_platform] } {
    puts "Seem to be not in a Unix-like environment."
    set ::my_platform "not_unix"
}
if {[string equal $tcl_platform(platform) "windows"] && \
	![string equal $::my_platform "Cygwin"] } {
    # we are on a plain Win32 system (without Cygwin)
    console show
}
if { [string equal $tcl_platform(platform) "unix"] } {
    # Something nicer than default, I think.
    ttk::style theme use clam
}
set ::debugFlag 1
# Get the directory of the script and assume that the
# executable programs are in the same place.
if { [string length [info script]] } {
    set ::binDir [file dirname [info script]]
} else {
    set ::binDir [file dirname [info nameofexecutable]]
}
set ::workDir [pwd]
cd $::workDir
if {$::debugFlag == 1} {
    puts "bin directory         : $::binDir"
    puts "current work directory: $::workDir"
}
set ::jobName ""
set ::scriptExt ".lua"
set ::subprocessId {}
set ::runProgramName [file join $::binDir e4shared]
set ::viewProgramName "paraview"

#-----------------------------------------------------------------------------
# Default options that may be reset by the user.

set optionsMain(maxWallClock) ""
set optionsMain(verbosity) 0
set optionsMain(noComplain) 0
set optionsMain(extras) ""

set optionsPost(tindxPlot) "all"
set optionsPost(outputFormat) "vtk-xml"
set optionsPost(addMach) 0
set optionsPost(addPitotP) 0
set optionsPost(addTotalP) 0
set optionsPost(extras) ""

set logFont {courier 12}

# -----------------------------------------------------------------------

proc startSubProcess { cmdString statusMessage } {
    if {$::debugFlag} { puts "cmdString is $cmdString" }
    appendToLogText "================ begin ================\n"
    set pipe [open "$cmdString"]
    fileevent $pipe readable [list Reader $pipe]
    set ::subprocessId [pid $pipe]
    updateStatusMessage $statusMessage
}; # end proc startSubProcess


proc Reader { pipe } {
    # Plumbing to enable this script to drive the external programs
    # which really do all of the work.
    if [eof $pipe] {
        catch {close $pipe}
	set ::subprocessId {}
	puts "End of pipe."
	flush stdout
	appendToLogText "================ end ==================\n"
	updateStatusMessage "Ready."
	return
    }
    gets $pipe line
    # puts $line; flush stdout
    appendToLogText "$line\n"; # reinstate the new-line character
}; # end proc Reader


proc killSubProcess {} {
    if { [string length $::subprocessId] > 0 } {
	appendToLogText "kill subprocess $::subprocessId\n"
	catch { exec kill $::subprocessId }
    }
}

# -----------------------------------------------------------------------

proc runPrep {} {
    set inputFile [file join $::workDir $::jobName.lua]
    if {[file exists $inputFile] == 0} {
	tk_messageBox -icon error -type ok \
	    -message "Could not find input script file: $inputFile"
	return
    }; # end if
    # Assume a unix environment.
    set cmd "|$::runProgramName --prep --job=$::jobName"
    startSubProcess $cmd "Running preparation stage, Please Wait..."
}; # end proc runPrep


proc helpPrep {} {
    # Assume a unix environment.
    set cmd "|$::runProgramName --help"
    startSubProcess $cmd "Usage message for preparation."
}; # end proc helpPrep


proc runMain {} {
    set cmd "|$::runProgramName --run --job=$::jobName"
    if { [string length [string trim $::optionsMain(maxWallClock)]] > 0 } {
	append cmd " --max-wall-clock=$::optionsMain(maxWallClock)"
    }
    set extras [string trim $::optionsMain(extras)]
    if { [string length $extras] > 0 } {
	append cmd " "
	append cmd $extras
    }
    startSubProcess $cmd "Running e4shared simulation, Please Wait..."
}; # end proc runMain


proc helpMain {} {
    set cmd "|$::runProgramName --help"
    startSubProcess $cmd "Usage message for e4shared."
}; # end proc helpMain


proc runPost {} {
    # Assume a unix environment.
    set cmd "|$::runProgramName --post --job=$::jobName"
    set addVarsStr ""
    if { $::optionsPost(addMach) } {
	append addVarsStr "mach,"
    }
    if { $::optionsPost(addPitotP) } {
	append addVarsStr "pitot,"
    }
    if { $::optionsPost(addTotalP) } {
	append addVarsStr "total-p,"
    }
    if [string length $addVarsStr] {
        append cmd " --add-vars=\"$addVarsStr\""
    }
    set tindxPlot [string trim $::optionsPost(tindxPlot)]
    if { ([string length $tindxPlot] > 0) } {
	# We are expecting one of "all", "last" or an integer+float pair
        # that has been pulled out of the job.times file.
	set indexValue [lindex [split $tindxPlot] 0]
	append cmd [join [list " --tindx-plot=" $indexValue] ""]
    }
    append cmd [join [list " --" $::optionsPost(outputFormat)] ""]
    set extras [string trim $::optionsPost(extras)]
    if { [string length $extras] > 0 } {
	append cmd " "
	append cmd $extras
    }
    startSubProcess $cmd "Running postprocessing stage, Please Wait..."
}; # end proc runPost


proc runViewer {} {
    set cmd "|$::viewProgramName"
    startSubProcess $cmd "Running Viewer, close that GUI window when finished..."
}; # end proc runViewer

# --------------------------------------------------------------

proc loadScriptFile { useExtEditor } {
    set fileName [tk_getOpenFile -title "Load script File..." \
        -initialdir $::workDir ]
    puts "loadScriptFile: fileName=$fileName"
    if {[string length $fileName] > 0} {
        set shortFileName [file tail $fileName]
	set ::scriptExt [file extension $shortFileName]
	set ::jobName [file root $shortFileName]
        set ::workDir [file dirname $fileName]; # save for next time
	if { $useExtEditor == 0 } {
	    readScriptFile $fileName
	} else {
	    set cmd "|/usr/bin/emacs $fileName"
	    set pipe [open "$cmd"]
	    fileevent $pipe readable [list Reader $pipe]
	    set ::editSubprocessId [pid $pipe]
	}
    }
}; # end proc loadScriptFile


proc readScriptFile { {fileName {}} } {
    if {[string length $fileName] == 0} {
	# Try to assemble a sensible filename from the available pieces.
	set fileName $::jobName
	append fileName $::scriptExt
	set fileName [file join $::workDir $fileName]
	puts "Assembled script file name: $fileName"
    }
    set fp [open $fileName "r"]
    $::scriptTextWidget delete 1.0 end
    $::scriptTextWidget insert end [read $fp]
    close $fp
}; # end proc readScriptFile


proc saveScriptFile {} {
    set fileName [tk_getSaveFile -title "Save script File..." \
        -initialdir $::workDir ]
    puts "saveScriptFile: fileName=$fileName"
    if {[string length $fileName] > 0} {
        set shortFileName [file tail $fileName]
	set ::scriptExt [file extension $shortFileName]
	set ::jobName [file root $shortFileName]
        set ::workDir [file dirname $fileName]; # save for next time
	set fp [open $fileName "w"]
	puts $fp [$::scriptTextWidget get 1.0 end-1c]
	close $fp
    }
}; # end proc loadScriptPyFile

# --------------------------------------------------------------

proc postMainMenu {} {
   # Set up menu bar for the main window
   # m0 is the top-level menubar
   # m1 is a first-level menu
   # m2 is a second-level menu
   set m0 [menu .menubar]
   . config -menu $m0
   wm title . "Eilmer4 Job Supervisor"

   set m1 [menu $m0.file -tearoff 0]
   $m0 add cascade -label File -menu $m1
   $m1 add command -label "Open script file (external editor)" \
       -command { loadScriptFile 1 }
   $m1 add command -label "Open script file (internal editor)" \
       -command { loadScriptFile 0; .menubar.file entryconfigure 2 -state normal }
   $m1 add command -label "Save script file (internal editor)" \
       -command { saveScriptFile } -state disabled
   $m1 add command -label "Clear log messages" \
       -command { clearLogText }
   $m1 add separator
   $m1 add command -label Quit -command "exit"

   set m1 [menu $m0.action -tearoff 1 -title "Action"]
   $m0 add cascade -label Action -menu $m1
   $m1 add command -label "Prepare grid and initial flow field" \
       -command {runPrep}
   $m1 add command -label "Run simulation" \
       -command {runMain}
   $m1 add command -label "Extract flow field data at specified time" \
       -command {runPost}
   $m1 add command -label "Run data-viewing program" \
       -command {runViewer}
   $m1 add separator
   $m1 add command -label "Kill subprocess" -command { killSubProcess }

   $m0 add command -label "Options..." \
       -command {displayOptionsDialog}

   set m1 [menu $m0.help -tearoff 0]
   $m0 add cascade -label Help -menu $m1
   $m1 add command -label "About..." \
       -command "showAboutBox"
   $m1 add command -label "e4shared usage" \
       -command { helpMain }
}; # end proc postMainMenu


proc initializeDisplay {} {
    postMainMenu
    # Items that make up the main window...
    #
    # (1) A line indicating the current jobName and extension.
    #
    set jobNameFrame [ttk::labelframe .jnf -text "Job indentity: " -labelanchor nw]
    pack [ttk::label $jobNameFrame.nl -text "Name:"] -side left
    set jobNameEntry [ttk::entry $jobNameFrame.jne -textvariable ::jobName -width 50]
    bind $jobNameEntry <Return> { readScriptFile }
    pack $jobNameEntry -side left -expand 1 -fill x
    pack [ttk::label $jobNameFrame.el -text "extension:"] -side left
    set extEntry [ttk::entry $jobNameFrame.ee -textvariable ::scriptExt -width 20 ]
    pack $extEntry -side left -expand 1 -fill x
    pack $jobNameFrame -side top -expand 1 -fill both -pady 5
    #
    # (2) A scrolling text window for the content of the user's script.
    #
    set scriptFrame [ttk::labelframe .tf0 -text "Input script:" -labelanchor nw]
    set ::scriptTextWidget [text $scriptFrame.t -height 20 -width 70 \
				-font $::logFont -wrap char \
		       -yscrollcommand [list $scriptFrame.vsb set] ]
    set scriptScrollBar [ttk::scrollbar $scriptFrame.vsb -orient vertical \
			     -command {$::scriptTextWidget yview} ]
    pack $::scriptTextWidget -side left -expand 1 -fill both
    pack $scriptScrollBar -side left -fill y
    pack $scriptFrame -side top -expand 1 -fill both -pady 5
    #
    # (3) A scrolling text window for the text output from running the programs.
    #
    set textFrame [ttk::labelframe .tf -text "Log of output:" -labelanchor nw]
    set ::logTextWidget [text $textFrame.t -height 10 -width 70 \
			     -font $::logFont -wrap char \
		       -yscrollcommand [list $textFrame.vsb set] ]
    set textScrollBar [ttk::scrollbar $textFrame.vsb -orient vertical \
			   -command {$::logTextWidget yview} ]
    pack $::logTextWidget -side left -expand 1 -fill both
    pack $textScrollBar -side left -fill y
    pack $textFrame -side top -expand 1 -fill both -pady 5
    #
    # (4) A status line giving some hint as to the current state.
    #
    set statusFrame [ttk::labelframe .fs -text "Status:" \
                         -labelanchor w -borderwidth 3 -relief flat]
    set statusEntry [ttk::entry $statusFrame.state -textvariable ::statusText ]
    pack $statusEntry -side right -expand 1 -fill x
    pack $statusFrame -side bottom -expand 1 -fill both -pady 5
}; # end proc initializeDisplay

proc clearLogText {} {
    $::logTextWidget delete 1.0 end
    update idletasks
}

proc appendToLogText { {newText ""} } {
    $::logTextWidget insert end $newText
    $::logTextWidget see end
    update idletasks
}

proc updateStatusMessage { {newText ""} } {
    set ::statusText $newText
    update idletasks
}

#-----------------------------------------------------------------------------

proc checkTimesFile {} {
    set fileName $::jobName
    append fileName ".times"
    set fp [open $fileName "r"]
    set tindxList {}
    foreach line [split [read -nonewline $fp] \n] {
	set tokens [split $line]
	if { [lindex $tokens 0] != "#" } {
	    lappend tindxList [join [list [lindex $tokens 0] [lindex $tokens 1]] " "]
	}
    }
    close $fp
    lappend tindxList "all"
    $::tindxListBox configure -values $tindxList
    return
}

proc displayOptionsDialog {} {
    if { [winfo exists .optionsDialog] } {
        # Dialog widget already exists; just show it
	# after 0.5 seconds to avoid other events bringing
	# the main window to the fore.
        after 500 {raise .optionsDialog}
    } else {
        # Create the dialog widget and show it.
        toplevel .optionsDialog
        set nb [ttk::notebook .optionsDialog.notebook]
        set bF [ttk::frame .optionsDialog.buttonFrame]
	#--------------------------------------------------
        $nb add [ttk::frame $nb.page1] -text "Prepare grids"
	set lab0 [ttk::label $nb.page1.lab0 -text "Use e4shared to generate grid and initial flow-field data."]
	set frame1 [ttk::frame $nb.page1.f1 -relief groove -borderwidth 2]
	set lab1 [ttk::label $frame1.lab1 -text "Standard options:"]
        set lab2 [ttk::label $frame1.lab2 -text "None, for the moment."]
	pack $lab1 $lab2 -anchor w
	pack $lab0 -anchor w -pady 10
	pack $frame1 -anchor w -fill x
	#--------------------------------------------------
        $nb add [ttk::frame $nb.page2] -text "Run simulation"
	set lab0 [ttk::label $nb.page2.lab0 -text "Use e4shared to generate subsequent flow-field data."]
	set fWallClock [ttk::frame $nb.page2.fWC]
	set labWallClock [ttk::label $fWallClock.labwc -text "Maximum Wall-Clock (sec):"]
	set entryWallClock [ttk::entry $fWallClock.entrywc -width 20 \
			     -textvariable ::optionsMain(maxWallClock)]
	pack $labWallClock -side left
	pack $entryWallClock -side left
	set labExtras [ttk::label $nb.page2.labe -text "Extra options:"]
	set entryExtras [ttk::entry $nb.page2.entrye -width 60 \
			     -textvariable ::optionsMain(extras)]
	pack $lab0 -anchor w -pady 10
        pack $fWallClock $labExtras $entryExtras -anchor w
	#--------------------------------------------------
        $nb add [ttk::frame $nb.page3] -text "Extract field data"
	set lab0 [ttk::label $nb.page3.lab0 -text "Use e4shared to extract flow-field data."]
	set frame1 [ttk::frame $nb.page3.f1]
	set labTime [ttk::label $frame1.labt -text "TimeIndex (integer, 'last' or 'all'):"]
	set ::tindxListBox [ttk::combobox $frame1.entryt -width 20 \
				-textvariable ::optionsPost(tindxPlot) \
				-values [list "all" "last" "0 0.0"] ]
	set checkTime [ttk::button $frame1.butt -text "Check times file" \
			   -command checkTimesFile]
	pack $labTime -side left
	pack $::tindxListBox -side left
	pack $checkTime -side left

	set buttonFrame [ttk::frame $nb.page3.bf]
	set formatFrame [ttk::labelframe $buttonFrame.ff -text "Output Format:" \
                             -labelanchor nw -borderwidth 2]
	set rb2 [ttk::radiobutton $formatFrame.rb2 \
		     -text "TECPLOT binary" \
		     -variable ::optionsPost(outputFormat) \
		     -value "tecplot"]
	set rb3 [ttk::radiobutton $formatFrame.rb3 \
		     -text "VTK-XML" \
		     -variable ::optionsPost(outputFormat) \
		     -value "vtk-xml"]
	set rb4 [ttk::radiobutton $formatFrame.rb4 \
		     -text "TECPLOT ascii" \
		     -variable ::optionsPost(outputFormat) \
		     -value "tecplot-ascii"]
	grid $rb2 -row 0 -column 0 -sticky w
	grid $rb3 -row 1 -column 0 -sticky w
	grid $rb4 -row 2 -column 0 -sticky w
	pack $formatFrame -side left -anchor w -padx 5 -pady 5

	set deFrame [ttk::labelframe $buttonFrame.lf -text "Add data elements:" \
			 -labelanchor nw -borderwidth 2]
	set cb1 [ttk::checkbutton $deFrame.cb1 \
		     -text "Mach number" \
		     -variable ::optionsPost(addMach)]
	set cb2 [ttk::checkbutton $deFrame.cb2 \
		     -text "Pitot pressure" \
		     -variable ::optionsPost(addPitotP)]
	set cb3 [ttk::checkbutton $deFrame.cb3 \
		     -text "Total pressure" \
		     -variable ::optionsPost(addTotalP)]
	grid $cb1 -row 0 -column 0 -sticky w
	grid $cb2 -row 1 -column 0 -sticky w
	grid $cb3 -row 2 -column 0 -sticky w
	pack $deFrame -side right -anchor e -padx 5 -pady 5

	set labExtras [ttk::label $nb.page3.labe -text "Extra options:"]
	set entryExtras [ttk::entry $nb.page3.entrye -width 60 \
			     -textvariable ::optionsPost(extras)]

	pack $lab0 -anchor w -pady 10 
        pack $frame1 $buttonFrame $labExtras $entryExtras -anchor w
	#--------------------------------------------------

        $nb select $nb.page1
        ttk::notebook::enableTraversal $nb

        set bDismiss [button $bF.dismiss -text "Dismiss or Lower this Window" \
            -command {lower .optionsDialog .} ]
        pack $bDismiss -side left -padx 4 -pady 4

        pack $nb $bF
    }; # end if
    wm title .optionsDialog "Options for actions"
}; # end proc

# --------------------------------------------------------

# Everything is ready; paint the main display...
initializeDisplay
displayOptionsDialog
lower .optionsDialog .
updateStatusMessage "Ready"

# --------------- end of e4console.tcl ---------------------
