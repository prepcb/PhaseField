  -  S   k820309    `          17.0        ��d                                                                                                           
       mpi_morton.F90 MPI_MORTON       D       NPTS_NEIGH PE_SOURCE MORTONBND R_MORTONBND MORTON_LIMITS IR_BUF IS_BUF NO_OF_MORTONBNDS_RECEIVED MORTON_LIMITS_SET MPI_TREE_SET COMMATRIX_SEND COMMATRIX_RECV COMMATRIX_GUARD COMMATRIX_PROL COMMATRIX_FLUX COMMATRIX_RESTRICT LADDRESS_GUARD LADDRESS_PROL LADDRESS_FLUX LADDRESS_RESTRICT EDGE_MARK NO_OF_DIAGONAL_EDGES TO_BE_SENT TO_BE_SENT_GUARD TO_BE_SENT_PROL TO_BE_SENT_FLUX TO_BE_SENT_RESTRICT TO_BE_RECEIVED LADD_STRT LADD_END TO_BE_RECEIVED_GUARD TO_BE_RECEIVED_PROL TO_BE_RECEIVED_FLUX TO_BE_RECEIVED_RESTRICT PE_SOURCE_GUARD PE_SOURCE_PROL PE_SOURCE_FLUX PE_SOURCE_RESTRICT MESSAGE_SIZE_CC MESSAGE_SIZE_FC MESSAGE_SIZE_EC MESSAGE_SIZE_NC MESSAGE_SIZE_WK MESS_SEGMENT_LOC TEMPRECV_BUF L_DATAPACKED LARGEST_NO_OF_BLOCKS LARGEST_NO_OF_BLOCKS_GUARD LARGEST_NO_OF_BLOCKS_PROL LARGEST_NO_OF_BLOCKS_FLUX LARGEST_NO_OF_BLOCKS_RESTRICT MAX_NO_TO_SEND MAX_NO_TO_SEND_GUARD MAX_NO_TO_SEND_PROL MAX_NO_TO_SEND_FLUX MAX_NO_TO_SEND_RESTRICT STRT_GUARD STRT_PROL STRT_FLUX STRT_RESTRICT NO_COMMATRIX_GUARD NO_COMMATRIX_PROL NO_COMMATRIX_FLUX NO_COMMATRIX_RESTRICT LPERIODICX LPERIODICY LPERIODICZ TREEINFO                                                     
                                                                                           �              3000                                                                             &                                                                                                                        &                   &                   &                                                                                                                        &                   &                   &                   &                                                                                                                        &                   &                   &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                       	                                                        
                                                                               @                                                                 &                                                      @                                                                 &                                                                                                                        &                                                                                                                        &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                                                        &                   &                                                                                            @                  p          p          p          p �          p          p          p �                                                                                                                                                             &                   &                   &                                                                                                                        &                   &                   &                                                                                                                        &                   &                   &                                                                                                                        &                   &                   &                                                                                                                         &                   &                   &                                                                                                                        &                   &                   &                                                                                                                         &                                                                                     !                                   &                                                                                     "                                   &                   &                   &                                                                                     #                                   &                   &                   &                                                                                     $                                   &                   &                   &                                                                                      %                                   &                   &                   &                                                                                     &                                   &                                                                                     '                                   &                                                                                     (                                   &                                                                                     )                                   &                                                                                       *     6                    p          p 6           p 6                                                                     +     6                    p          p 6           p 6                                                                     ,     6                    p          p 6           p 6                                                                     -     6                    p          p 6           p 6                                                                     .     6                    p          p 6           p 6                                                                   /                                   &                                                                                     0                   
                &                                                                                       1                         p          p            p                                                                      2                                                        3                                                        4                                                        5                                                        6                                                        7                                                        8                                                        9                                                        :                                                        ;                                                        <                                                        =                                                        >                                                        ?                                                        @                                                        A                                                        B                                                        C                                                        D                                                        E                                                        F                              @                           G     '�              
      #COORD H   #BSIZE I   #BND_BOX J   #PARENT K   #WHICH_CHILD L   #NEWCHILD M   #NEIGH N   #LREFINE O   #NODETYPE P   #EMPTY Q                � $                              H                              
  p          p            p                                       � $                              I                             
  p          p            p                                       � $                              J            0                 
  p          p          p            p          p                                       � $                              K            `                   p          p            p                                       � $                              L     h                          � $                              M     l                          � $                              N            p                   p          p          p            p          p                                       � $                              O     �                          � $                              P     �       	                   � $                              Q     �       
         �   "      fn#fn     �   X  b   uapp(MPI_MORTON $     @   J  PARAMESH_DIMENSIONS    Z  t       NPTS_NEIGH    �  �       PE_SOURCE    Z  �       MORTONBND      �       R_MORTONBND    �  �       MORTON_LIMITS    �  �       IR_BUF    b	  �       IS_BUF *   
  @       NO_OF_MORTONBNDS_RECEIVED "   F
  @       MORTON_LIMITS_SET    �
  @       MPI_TREE_SET    �
  �       PE_REMOTE    R  �       PE_DESTINATION    �  �       COMMATRIX_SEND    j  �       COMMATRIX_RECV     �  �       COMMATRIX_GUARD    �  �       COMMATRIX_PROL    >  �       COMMATRIX_FLUX #   �  �       COMMATRIX_RESTRICT    �  �       LADDRESS_GUARD    *  �       LADDRESS_PROL    �  �       LADDRESS_FLUX "   r  �       LADDRESS_RESTRICT      �       EDGE_MARK %   �  @       NO_OF_DIAGONAL_EDGES    *  �       TO_BE_SENT !   �  �       TO_BE_SENT_GUARD     �  �       TO_BE_SENT_PROL     ^  �       TO_BE_SENT_FLUX $     �       TO_BE_SENT_RESTRICT    �  �       TO_BE_RECEIVED    �  �       LADD_STRT      �       LADD_END %   �  �       TO_BE_RECEIVED_GUARD $   f  �       TO_BE_RECEIVED_PROL $   "  �       TO_BE_RECEIVED_FLUX (   �  �       TO_BE_RECEIVED_RESTRICT     �  �       PE_SOURCE_GUARD    &  �       PE_SOURCE_PROL    �  �       PE_SOURCE_FLUX #   >  �       PE_SOURCE_RESTRICT     �  �       MESSAGE_SIZE_CC     ^  �       MESSAGE_SIZE_FC     �  �       MESSAGE_SIZE_EC     �  �       MESSAGE_SIZE_NC        �       MESSAGE_SIZE_WK !   �   �       MESS_SEGMENT_LOC    :!  �       TEMPRECV_BUF    �!  �       L_DATAPACKED %   Z"  @       LARGEST_NO_OF_BLOCKS +   �"  @       LARGEST_NO_OF_BLOCKS_GUARD *   �"  @       LARGEST_NO_OF_BLOCKS_PROL *   #  @       LARGEST_NO_OF_BLOCKS_FLUX .   Z#  @       LARGEST_NO_OF_BLOCKS_RESTRICT    �#  @       MAX_NO_TO_SEND %   �#  @       MAX_NO_TO_SEND_GUARD $   $  @       MAX_NO_TO_SEND_PROL $   Z$  @       MAX_NO_TO_SEND_FLUX (   �$  @       MAX_NO_TO_SEND_RESTRICT    �$  @       STRT_GUARD    %  @       STRT_PROL    Z%  @       STRT_FLUX    �%  @       STRT_RESTRICT #   �%  @       NO_COMMATRIX_GUARD "   &  @       NO_COMMATRIX_PROL "   Z&  @       NO_COMMATRIX_FLUX &   �&  @       NO_COMMATRIX_RESTRICT    �&  @       LPERIODICX    '  @       LPERIODICY    Z'  @       LPERIODICZ    �'  �       TREEINFO    i(  �   a   TREEINFO%COORD    )  �   a   TREEINFO%BSIZE !   �)  �   a   TREEINFO%BND_BOX     ]*  �   a   TREEINFO%PARENT %   �*  H   a   TREEINFO%WHICH_CHILD "   A+  H   a   TREEINFO%NEWCHILD    �+  �   a   TREEINFO%NEIGH !   E,  H   a   TREEINFO%LREFINE "   �,  H   a   TREEINFO%NODETYPE    �,  H   a   TREEINFO%EMPTY 