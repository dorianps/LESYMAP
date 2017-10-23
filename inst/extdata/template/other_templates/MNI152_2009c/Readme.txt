ICBM 2009c Nonlinear Symmetric

Obtained in .nii format from:
http://nist.mni.mcgill.ca/?p=904

Files were then opened in ITKsnap and saved in .nii.gz format to save space.

The _skullnoface.nii.gz file contains a mask with that includes the skull but
does not include the face. It is used to improve registrations when registering 
without skull stripping, and is used during skull-stipping itself to improve the result.
This mask was created by registering another existing template on the ICBM 2009c
and bringing the mask in ICBM space.
