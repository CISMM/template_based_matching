#############################################################################
# Sets up the control panel and functions for the Template-based matching program.

# Global variable to remember where they are saving files.
set fileinfo(open_dir) ""
#set fileinfo(open_dir) "C:\\"

# Global variable for background color
set background_color #000000
set background_r 0
set background_g 0
set background_b 0

###########################################################
# Add a button to let the user quit the program.

button .quitbutton -text "Quit" -command {set quit 1}
pack .quitbutton -fill x

###########################################################
# Add a button to let the user add a new data set to the
# ones that have been opened already.  Don't quit the
# program if the user presses cancel.

#button .newdata -text "Change DataSet" -command {ask_user_for_data_filename; set quit 0}
#pack .newdata -fill x

###########################################################
# Add a button to let the user add a new template to the
# ones that have been opened already.  Don't quit the
# program if the user presses cancel.

button .newtemplate -text "Add Template" -command {ask_user_for_template_filename; set quit 0}
pack .newtemplate -fill x

###########################################################
# Add a min/max slider to threshold the volumes in the
# component clusters.

label .volumethreshlabel -text "Volume Thresholds"
pack .volumethreshlabel
minmaxscale .volumethresh 0 1 50 min_norm_volume max_norm_volume
pack .volumethresh
set min_norm_volume 0
set max_norm_volume 1

###########################################################
# Add a button to let the user optimize and save.
# Put a trace in place so that it calls the save dialog to
# pick a save file name when it is actived.

button .savebutton -text "Save" -command {set save 1}
pack .savebutton -fill x -side bottom
trace variable save w save_changed

proc save_changed { varName index op } {
    global save savefilename fileinfo

    if {$save == 1} {
	set types { {"Comma-Separated Value Reports" "*.csv"} }
	set filename [tk_getSaveFile -filetypes $types \
		-initialdir $fileinfo(open_dir) \
		-defaultextension ".csv" \
		-title "Name for statistics file"]
	if {$filename != ""} {
	    # setting this variable triggers a callback in C code
	    # which opens the file.
	    # dialog check whether file exists.
	    set savefilename $filename
	    set fileinfo(open_dir) [file dirname $filename]
	}
    } else {
	set logfilename ""
    }
}

###########################################################
# Ask user for the name of the data image file they want to open,
# or else set the quit value.  The variable to set for the
# name is "data_filename".

set save_template_filename ""
proc ask_user_for_save_template_filename { } {
	global save_template_filename quit fileinfo
		
	set types { {"UNCA files" "*.unca"} }
	set save_template_filename [tk_getSaveFile -filetypes $types \
		-defaultextension ".unca" \
		-initialdir $fileinfo(open_dir) \
		-title "Specify a filename to save template to"]
	# If we don't have a name, quit.
	if {$save_template_filename == ""} {
		set quit 1
	} else {
	  # Look in the same directory for files next time
        set fileinfo(open_dir) [file dirname $save_template_filename]
	}
}

###########################################################
# Pop up a text window displaying progress in the matching.
# Include a "Cancel" button so that people can break out in
# the middle.

proc track_progress { } {
	global comp_id comp_count temp_id temp_count save_cancel
#	global angle

	toplevel .progress

	frame .progress.comp
	pack .progress.comp

	label .progress.comp.component -text "Testing component "
	pack .progress.comp.component -side left
	label .progress.comp.id -textvariable comp_id
	pack .progress.comp.id -side left
	label .progress.comp.of -text " of "
	pack .progress.comp.of -side left
	label .progress.comp.total -textvariable comp_count
	pack .progress.comp.total -side left

	frame .progress.temp
	pack .progress.temp

	label .progress.temp.template -text "Testing template "
	pack .progress.temp.template -side left
	label .progress.temp.id -textvariable temp_id
	pack .progress.temp.id -side left
	label .progress.temp.of -text " of "
	pack .progress.temp.of -side left
	label .progress.temp.total -textvariable temp_count
	pack .progress.temp.total -side left

#	frame .progress.angle
#	pack .progress.angle

#	label .progress.angle.template -text "Testing angle "
#	pack .progress.angle.template -side left
#	label .progress.angle.id -textvariable angle
#	label .progress.angle.of -text " of 360"
#	pack .progress.angle.of -side left

	set save_cancel 0
	button .progress.cancelbutton -text "Cancel" -command {set save_cancel 1}
	pack .progress.cancelbutton -fill x
}

proc done_tracking_progress { } {
	destroy .progress
}

###########################################################
# Add a background-color selector button along with a swatch
# of color that indicates the selected color.
#frame .backcolor
#button .backcolor.set_color -text "Set background color" -command {
#    choose_color background_color "Choose background color" .backcolor
#    .backcolor.colorsample configure -bg $background_color
#    set_background_color
#}
#
# This sample frame displays the color of the background
#button .backcolor.colorsample -relief groove -bd 2 -bg $background_color \
#        -command { .backcolor.set_color invoke}
#
#pack .backcolor.set_color -side left
#pack .backcolor.colorsample
#pack .backcolor
#
#proc set_background_color {} {
#    global background_changed background_color background_r background_g #background_b
#    # Extract three component colors of background_color 
#    # and save into background_r g b
#    scan $background_color #%02x%02x%02x background_r background_g #background_b
#    set background_changed 1
#}
#

###########################################################
# Ask user for the name of the data image file they want to open,
# or else set the quit value.  The variable to set for the
# name is "data_filename".

set data_filename ""
proc ask_user_for_data_filename { } {
	global data_filename quit fileinfo
		
	set types { {"All Files" "*.*"} }
	set data_filename [tk_getOpenFile -filetypes $types \
		-defaultextension ".tif" \
		-initialdir $fileinfo(open_dir) \
		-title "Specify a data file to display"]
	# If we don't have a name, quit.
	if {$data_filename == ""} {
		set quit 1
	} else {
	  # Look in the same directory for files next time
        set fileinfo(open_dir) [file dirname $data_filename]
	}
}

###########################################################
# Ask user for the name of the template image file they want to open,
# or else set the quit value.  The variable to set for the
# name is "template_filename".

set template_filename ""
proc ask_user_for_template_filename { } {
	global template_filename quit fileinfo
		
	set types { {"All Files" "*.*"} }
	set template_filename [tk_getOpenFile -filetypes $types \
		-defaultextension ".tif" \
		-initialdir $fileinfo(open_dir) \
		-title "Specify a template file to display"]
	# If we don't have a name, quit.
	if {$template_filename == ""} {
		set quit 1
	} else {
	  # Look in the same directory for files next time
        set fileinfo(open_dir) [file dirname $template_filename]
	}
}

###########################################################
# Creates a standard dialog for choosing a color, and sets the
# value of a global variable based on the user's choice. If the
# user hits "cancel", the global variable's value is unchanged. 
# Returns 1 if user hit OK, and 0 if user hit Cancel
proc choose_color { color_var_name {title "Choose color"} {parent .} } {
    global nmInfo
    upvar #0 $color_var_name color_var
    set color [tk_chooseColor -title "$title" \
	    -initialcolor $color_var -parent $parent]
    if { $color != "" } {
	set color_var $color
        return 1
    }
    return 0
}