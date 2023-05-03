# Percy++ Copyright 2007,2012,2013,2014
# Ian Goldberg <iang@cs.uwaterloo.ca>,
# Casey Devet <cjdevet@uwaterloo.ca>
##
# This program is free software; you can redistribute it and/or modify
# it under the terms of version 2 of the GNU General Public License as
# published by the Free Software Foundation.
##
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
##
# There is a copy of the GNU General Public License in the COPYING file
# packaged with this plugin; if you cannot find it, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301 USA

import paramiko
import os
import sys

devnull = open(os.devnull, 'wb')


class SSHConfig (object):
    def __init__(self, files=[]):
        self.hosts = paramiko.SSHConfig()
        for f in files:
            try:
                with open(os.path.expanduser(f), 'r') as config:
                    self.hosts.parse(config)
            except:
                print >>sys.stderr, 'Warning: could not read config file:', f

    def lookup(self, host):
        return self.hosts.lookup(host)


class RemoteCommand (object):
    def __init__(self, host, command, directory=None):
        self._host = host
        cmd = 'echo $$; exec %s' % (command)
        if directory != None:
            cmd = 'cd %s; ' % (directory) + cmd
        self._chan = host.client._transport.open_session()
        self._chan.exec_command(cmd)
        self.stdin = self._chan.makefile('wb')
        self.stdout = self._chan.makefile('rb')
        self.stderr = self._chan.makefile_stderr('rb')
        self.pid = int(self.stdout.readline().strip())

    def is_running(self):
        return not self._chan.exit_status_ready()

    def kill(self):
        if self.is_running():
            self._host.client.exec_command('kill -- -%d' % (self.pid))

    def wait(self):
        return self._chan.recv_exit_status()

    def get_exit_code(self):
        return self._chan.exit_status


class RemoteHost (object):
    _timeout = 5

    def __init__(self, host, config={}):
        self.host = host
        self.config = config
        self.addr = config.get('hostname', host)
        self.user = config.get('user', None)
        self.client = paramiko.SSHClient()
        self.client.load_system_host_keys()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.client.connect(self.addr,
                            port=config.get('port', 22),
                            username=self.user,
                            password=config.get('password', None),
                            key_filename=config.get('identifyfile', None),
                            timeout=self._timeout)

    def exec_command(self, command, directory=None):
        if isinstance(command, str):
            cmd = command
        elif isinstance(command, list):
            cmd = ' '.join(
                map(lambda x: '"' + x + '"' if (' ' in x or x == '') else x, command))
        else:
            raise TypeError('command must be a str or list.')
        return RemoteCommand(self, cmd, directory)

    def close(self):
        self.client.close()

    def __str__(self):
        ret = self.host
        if self.user != None:
            ret = self.user + '@' + ret
        return ret
