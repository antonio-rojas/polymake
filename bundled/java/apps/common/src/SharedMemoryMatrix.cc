/* Copyright (c) 1997-2018
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/common/SharedMemoryMatrix.h"
#include <sys/ipc.h>
#include <sys/shm.h>
#include <string>
#include <cerrno>

namespace polymake { namespace common {

void SharedMemorySegment::resize(size_t size)
{
   shmid=shmget(IPC_PRIVATE, size, 0600);
   if (shmid<0)
      throw std::runtime_error("shmget error " + std::to_string(errno));

   void* p=shmat(shmid, NULL, 0);
   if (p==(void*)-1L) {
      std::string err="shmat error " + std::to_string(errno);
      shmctl(shmid, IPC_RMID, 0);
      throw std::runtime_error(err);
   }
   shmaddr.set(p);
}

SharedMemorySegment::~SharedMemorySegment()
{
   if (void* p=shmaddr.get()) {
      shmdt(p);
      shmctl(shmid, IPC_RMID, 0);
   }
}

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
