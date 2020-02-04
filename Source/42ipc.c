/*    This file is distributed with 42,                               */
/*    the (mostly harmless) spacecraft dynamics simulation            */
/*    created by Eric Stoneking of NASA Goddard Space Flight Center   */

/*    Copyright 2010 United States Government                         */
/*    as represented by the Administrator                             */
/*    of the National Aeronautics and Space Administration.           */

/*    No copyright is claimed in the United States                    */
/*    under Title 17, U.S. Code.                                      */

/*    All Other Rights Reserved.                                      */


#include "42.h"
#include <time.h>
#include <unistd.h>

#ifdef _ENABLE_GMSEC_
   #include "gmseckit.h"
   GMSEC_Config cfg;
   GMSEC_ConnectionMgr ConnMgr;
   GMSEC_Status status;
   void WriteToGmsec(GMSEC_ConnectionMgr ConnMgr,GMSEC_Status status);
   void ReadFromGmsec(GMSEC_ConnectionMgr ConnMgr,GMSEC_Status status);
#endif

/* #ifdef __cplusplus
** namespace _42 {
** using namespace Kit;
** #endif
*/

void WriteToFile(FILE *StateFile, char **Prefix, long Nprefix, long EchoEnabled);
void WriteToSocket(SOCKET Socket,  char **Prefix, long Nprefix, long EchoEnabled);
void ReadFromFile(FILE *StateFile, long EchoEnabled);
void ReadFromSocket(SOCKET Socket, long EchoEnabled);

/*********************************************************************/
void InitInterProcessComm(void)
{
      FILE *infile;
      char junk[120],newline,response[120];
      struct IpcType *Ipc;
      long Iipc,Ipx;
      char FileName[80],Prefix[80];

      infile = FileOpen(InOutPath,"Inp_IPC.txt","rt");
      fscanf(infile,"%[^\n] %[\n]",junk,&newline);
      fscanf(infile,"%ld %[^\n] %[\n]",&Nipc,junk,&newline);
      IPC = (struct IpcType *) calloc(Nipc,sizeof(struct IpcType));
      for(Iipc=0;Iipc<Nipc;Iipc++) {
         Ipc = &IPC[Iipc];
         fscanf(infile,"%[^\n] %[\n]",junk,&newline);
         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
         Ipc->Mode = DecodeString(response);
         fscanf(infile,"%ld %[^\n] %[\n]",&Ipc->AcsID,junk,&newline);
         fscanf(infile,"\"%[^\"]\" %[^\n] %[\n]",FileName,junk,&newline);
         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
         Ipc->SocketRole = DecodeString(response);
         fscanf(infile,"%s %ld %[^\n] %[\n]",Ipc->HostName,&Ipc->Port,junk,&newline);
         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
         Ipc->AllowBlocking = DecodeString(response);
         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
         Ipc->EchoEnabled = DecodeString(response);
         fscanf(infile,"%ld %[^\n] %[\n]",&Ipc->Nprefix,junk,&newline);
         Ipc->Prefix = (char **) calloc(Ipc->Nprefix,sizeof(char *));
         for(Ipx=0;Ipx<Ipc->Nprefix;Ipx++) {
            fscanf(infile,"\"%[^\"]\" %[^\n] %[\n]",Prefix,junk,&newline);
            Ipc->Prefix[Ipx] = (char *) calloc(strlen(Prefix)+1,sizeof(char));
            strcpy(Ipc->Prefix[Ipx],Prefix);
         }
         
         Ipc->Init = 1;

         if (Ipc->Mode == IPC_TX) {
            if (Ipc->SocketRole == IPC_SERVER) {
               Ipc->Socket = InitSocketServer(Ipc->Port,Ipc->AllowBlocking);
            }
            else if (Ipc->SocketRole == IPC_CLIENT) {
               Ipc->Socket = InitSocketClient(Ipc->HostName,Ipc->Port,Ipc->AllowBlocking);
            }
            #ifdef _ENABLE_GMSEC_
            else if (I->SocketRole == IPC_GMSEC_CLIENT) {
               status = statusCreate();
               cfg = configCreate();
               ConnMgr = ConnectToMBServer(I->HostName,I->Port,status,cfg);
               connectionManagerSubscribe(ConnMgr,"GMSEC.42.RX.>",status);
               CheckGmsecStatus(status);
            }
            #endif
            else {
               printf("Oops.  Unknown SocketRole %ld for IPC[%ld] in InitInterProcessComm.  Bailing out.\n",Ipc->SocketRole,Iipc);
               exit(1);
            }
         }
         else if (Ipc->Mode == IPC_RX) {
            if (Ipc->SocketRole == IPC_SERVER) {
               Ipc->Socket = InitSocketServer(Ipc->Port,Ipc->AllowBlocking);
            }
            else if (Ipc->SocketRole == IPC_CLIENT) {
               Ipc->Socket = InitSocketClient(Ipc->HostName,Ipc->Port,Ipc->AllowBlocking);
            }
            #ifdef _ENABLE_GMSEC_
            else if (I->SocketRole == IPC_GMSEC_CLIENT) {
               status = statusCreate();
               cfg = configCreate();
               ConnMgr = ConnectToMBServer(I->HostName,I->Port,status,cfg);
               connectionManagerSubscribe(ConnMgr,"GMSEC.42.TX.>",status);
               CheckGmsecStatus(status);
            }
            #endif
            else {
               printf("Oops.  Unknown SocketRole %ld for IPC[%ld] in InitInterProcessComm.  Bailing out.\n",Ipc->SocketRole,Iipc);
               exit(1);
            }

         }
         else if (Ipc->Mode == IPC_ACS) {
            Ipc->Socket = InitSocketServer(Ipc->Port,Ipc->AllowBlocking);
         }
         else if (Ipc->Mode == IPC_WRITEFILE) {
            Ipc->File = FileOpen(InOutPath,FileName,"wt");
         }
         else if (Ipc->Mode == IPC_READFILE) {
            Ipc->File = FileOpen(InOutPath,FileName,"rt");
         }
         else if (Ipc->Mode == IPC_FFTB) {
            Ipc->SocketRole = IPC_CLIENT; /* Spirent is Host */
            Ipc->Socket = InitSocketClient(Ipc->HostName,Ipc->Port,Ipc->AllowBlocking);
         }
      }
      fclose(infile);
}
/*********************************************************************/
void InterProcessComm(void)
{
      struct IpcType *Ipc;
      long Iipc;
      
      for(Iipc=0;Iipc<Nipc;Iipc++) {
         Ipc = &IPC[Iipc];
         if (Ipc->Mode == IPC_TX) {
            if (Ipc->SocketRole != IPC_GMSEC_CLIENT) {
               WriteToSocket(Ipc->Socket,Ipc->Prefix,Ipc->Nprefix,Ipc->EchoEnabled);
            }
            #ifdef _ENABLE_GMSEC_
            else {
               WriteToGmsec(ConnMgr,status,I->Prefix,I->Nprefix,I->EchoEnabled);
            }
            #endif
         }
         else if (Ipc->Mode == IPC_RX) {
            if (Ipc->SocketRole != IPC_GMSEC_CLIENT) {
               ReadFromSocket(Ipc->Socket,Ipc->EchoEnabled);
            }
            #ifdef _ENABLE_GMSEC_
            else {
               ReadFromGmsec(ConnMgr,status,I->EchoEnabled);
            }
            #endif
         }
         else if (Ipc->Mode == IPC_WRITEFILE) {
            WriteToFile(Ipc->File,Ipc->Prefix,Ipc->Nprefix,Ipc->EchoEnabled);
         }
         else if (Ipc->Mode == IPC_READFILE) {
            ReadFromFile(Ipc->File,Ipc->EchoEnabled);
         }
         #ifdef _ENABLE_FFTB_CODE_
         else if (I->Mode == IPC_FFTB) {
            SendStatesToSpirent();
         }
         #endif
      }
}

/* #ifdef __cplusplus
** }
** #endif
*/
