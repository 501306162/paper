function Dd = D_new(dr,dc,dz)  %2-D
if(nargin<3)  %2d

      [dr_dc,dr_dr]=gradient(dr);
      [dc_dc,dc_dr]=gradient(dc);
      
      Dd=[dr_dr(:)';dr_dc(:)';
          dc_dr(:)';dc_dc(:)'];
else    %3d
      [dr_dc,dr_dr,dr_dz]=gradient(dr);
      [dc_dc,dc_dr,dc_dz]=gradient(dc);
      [dz_dc,dz_dr,dz_dz]=gradient(dz);
      Dd=[dr_dr(:)';dr_dc(:)';dr_dz(:)';
          dc_dr(:)';dc_dc(:)';dc_dz(:)';
          dz_dr(:)';dz_dc(:)';dz_dz(:)';];
end
end