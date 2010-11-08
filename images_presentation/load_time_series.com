gfx read nodes Reentry.exnode
gfx read elem Reentry.exelem generate_faces_and_lines

$number_of_frames = 100;
for ($file_number = 0; $file_number<=$number_of_frames; $file_number++)
{
   gfx read nodes Reentry_$file_number.exnode time $file_number
}
