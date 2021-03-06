#  Copyright (c) 1997-2018
#  Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
#  http://www.polymake.org
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation; either version 2, or (at your option) any
#  later version: http://www.gnu.org/licenses/gpl.txt.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#-------------------------------------------------------------------------------

use strict;
use namespaces;
use warnings qw(FATAL void syntax misc);
use feature 'state';

package Polymake::User;

declare $application;

sub application {
   if (@_>1) {
      die "usage: application [ \"name\" ]\n";
   } elsif (my ($new_app)=@_) {

      # This magic provides automatic loading of applications when they are first mentioned
      # as a prefix of a user function in the shell input, in a script, documentation example, or tutorial.
      state $register_autoload=namespaces::set_autolookup(\&Core::Application::try_auto_load);

      if (defined wantarray) {
         if (ref($new_app)) {
            warn_print( "application() call without effect as the application ", $new_app->name, " already loaded" );
            $new_app;
         } else {
            add Core::Application($new_app);
         }
      } else {
         $new_app=add Core::Application($new_app) unless ref($new_app);
         if (defined($application)) {
            return if $application == $new_app;
            readwrite($application);
         }
         $application=$new_app;
         readonly($application);
      }
   } else {
      # tell the current application
      $application;
   }
}

#################################################################################
sub include {
   my $rc=$application->include_rules(@_);
   if ($rc != @_) {
      foreach (@_) {
         my ($app, $rulefile)= /^($id_re)::/o ? ($application->used->{$1}, $') : ($application, $_);
         my ($filename, $ext, $rule_key, $rc)=$app->lookup_rulefile($rulefile);
         if (!$rc) {
            if ($app->configured->{$rule_key} < 0) {
               warn_print( "rulefile $_ is disabled by auto-configuration.\n",
                           "Try  reconfigure \"$_\";  if you really need it." );
            } elsif ($app->configured->{$rule_key} =~ /^0\#(?=.)/) {
               warn_print( "rulefile $_ can't be loaded because of unsatisfied dependency on $'\n",
                           "Try to reconfigure the prerequisites if you really need it." );
            }
         }
      }
   }
}
#################################################################################
sub load {
   my ($name)=@_;
   my $filename=$name;
   replace_special_paths($filename);
   unless (-f $filename
           or
           defined($application->default_file_suffix)  &&  $filename !~ /\.\w+$/  &&
           -f ($filename .= "." . $application->default_file_suffix)) {
      croak( "no such file: $name" );
   }
   my $obj=load Core::Object($filename);
   $obj;
}
#################################################################################
sub save {
   my ($obj, $filename)=@_;
   if (!is_object($obj) || !instanceof Core::Object($obj)) {
      croak( "only polymake Objects can be stored with the function save()" );
   }
   if (defined($obj->parent)) {
      croak( "A sub-object can't be saved without its parent" );
   }
   if (!$obj->persistent || defined($filename)) {
      if (defined $filename) {
         replace_special_paths($filename);
         if ($filename !~ /\.\w+$/) {
            $filename .= "." . $obj->default_file_suffix;
         }

      } else {
         if (!length($obj->name)) {
            $obj->name=Core::name_of_arg_var(0);
         }
         if (length($obj->name)) {
            $filename=$obj->name.".".$obj->default_file_suffix;
            if (-f $filename) {
               if ($Shell->interactive) {
                  print "The file $filename already exists.\n",
                        "Please specify another file name or confirm it to be overwritten:\n";
                  $filename=$Shell->enter_filename($filename, { prompt => "data file" }) or return;
               } else {
                  croak( "Can't save an object '", $obj->name, "' since the file $filename already exists.\n",
                         "Please specify the explicit file name as the second argument to save() or delete the existing file (unlink \"$filename\")" );
               }
            }
         } else {
            if ($Shell->interactive) {
               print "Please specify the file name for the anonymous ", $obj->type->full_name, " object:\n";
               $filename=$Shell->enter_filename("", { prompt => "data file" }) or return;
            } else {
               my $i=1;
               while (-f ($filename=($obj->name=$obj->type->name."_$i".".".$obj->default_file_suffix))) { ++$i }
               warn_print( "saving object as $filename" );
            }
         }
      }
      $obj->persistent=new Core::XMLfile($filename);
      $obj->changed=1;

   } elsif (! $obj->changed) {
      warn_print( "no changes need to be saved" );
      return;
   }

   $obj->save;
}
#################################################################################
sub save_schema {
   if (@_<2) {
      die "usage: save_schema(Object or ObjectType, ..., \"filename\"\n";
   }
   my $filename=pop;
   if ($filename !~ /\.\w+$/) {
      $filename .= ".rng";
   }
   replace_special_paths($filename);
   my $xf=new Core::XMLfile($filename);
   $xf->save_schema(map {
           if (is_object($_) && UNIVERSAL::can($_, "type")) {
              $_->type
           } else {
              die "argument ", ref($_) || "'$_'", " is not an object or a type expression\n";
           }
        } @_);
}
#################################################################################
sub load_data {
   my ($filename)=@_;
   replace_special_paths($filename);
   my $xf=new Core::XMLfile($filename);
   scalar($xf->load_data);
}

sub save_data {
   my ($data, $filename, $description)=@_;
   if (!is_object($data)) {
      croak( "only complex data types can be saved in separate files" );
   }
   if (instanceof Core::Object($data)) {
      croak( "an object of type ", $data->type->full_name, " can't be saved with save_data: please use save() instead" );
   }
   if (defined (my $type=UNIVERSAL::can($data, ".type"))) {
      $type=$type->();
      if (instanceof Core::PropertyType($type)) {
         replace_special_paths($filename);
         my $xf=new Core::XMLfile($filename);
         $xf->save_data($data, $description);
         return;
      }
   }
   croak( "can't save ", ref($data), ": it does not belong to any declared property type" );
}
#################################################################################
use subs qw(rename unlink mkdir chdir rmdir);

sub rename { my ($from, $to)=@_; replace_special_paths($from, $to); CORE::rename($from, $to) or die "rename failed: $!\n"; }
sub unlink { my @list=@_; replace_special_paths(@list); CORE::unlink(@list) or die "unlink failed: $!\n"; }
sub mkdir { my ($path, $mask)=@_; replace_special_paths($path); CORE::mkdir($path, $mask || 0755) or die "mkdir failed: $!\n"; }
sub rmdir { my ($path)=@_; replace_special_paths($path); CORE::rmdir($path) or die "rmdir failed: $!\n"; }

sub chdir {
   my $path;
   if (my ($path)=@_) {
      replace_special_paths($path) if is_string($path);
      CORE::chdir($path) or die "chdir failed: $!\n";
   } else {
      CORE::chdir;
   }
}

sub pwd { print Cwd::cwd }

#################################################################################
sub prefer {
   if ($_[0] =~ /^($id_re)::/o) {
      shift;  application($1)->prefer($', @_);
   } else {
      $application->prefer(@_);
   }
}

sub prefer_now {
   if ($_[0] =~ /^($id_re)::/o) {
      shift;  application($1)->prefer_now($', @_);
   } else {
      $application->prefer_now(@_);
   }
}

# an alias, for the sake of symmetry
*set_preference=\&prefer;

sub reset_preference {
   if ($_[0] =~ /^($id_re)::/o) {
      shift;  application($1)->reset_preference($', @_);
   } else {
      $application->reset_preference(@_);
   }
}
#################################################################################
sub disable_rules {
   $application->disable_rules(@_);
}
#################################################################################
sub set_custom {
   $application->_set_custom($Core::Prefs->custom, Core::name_of_custom_var(1));
}

sub reset_custom {
   $application->_reset_custom($Core::Prefs->custom, Core::name_of_custom_var(0));
}
#################################################################################
sub script {
   my $name=shift;
   replace_special_paths($name);
   my ($code, $full_path, $in_app)=find Core::StoredScript($name);
   if (defined $code) {
      local *ARGV=\@_;
      &$code;
   } else {
      local @ARGV = @_;
      local if (defined($in_app)) {
         if ($in_app != $User::application) {
            local $User::application = $in_app;
            if (ref($INC[0])) {
               local $INC[0] = $in_app;
            } else {
               local unshift @INC, $in_app;
            }
         }
      } elsif (defined (local scalar $User::application)) {
         if (!ref($INC[0])) {
            local unshift @INC, $User::application;
         }
      } else {
         if (!ref($INC[0])) {
            local unshift @INC, new Core::NeutralScriptLoader();
         }
      }
      $name="script" . (defined($in_app) ? ":" : "=") . $full_path;
      if (wantarray) {
         my @ret=do $name;
         die $@ if $@;
         @ret
      } elsif (defined wantarray) {
         my $ret=do $name;
         die $@ if $@;
         $ret
      } else {
         do $name;
         die $@ if $@;
      }
   }
}
#################################################################################
# print boolean values in legible form: true and false instead of 1 and empty string
# enforce creation of a unique lexical scope with this operation inherited by all nested packages
use namespaces 'Polymake::User';
namespaces::memorize_lexical_scope;
namespaces::intercept_operation(undef, "P", "bool");

#################################################################################
# prepare for custom variables and preferences

package Polymake::User::Verbose;
*Polymake::Verbose::=get_symtab(__PACKAGE__);

Core::add_custom_vars sub {
   my $ch = $Core::Prefs->create_custom("Polymake::User");

   $ch->pkg_help->{__PACKAGE__}=<<'.';
# The following variables control the display of various informational message classes.
.
   declare $credits=1;
   $ch->add('$credits', <<'.');
# Display the copyright notices and URLs of third-party software:
# 0 - never, 1 - at the first use in the current session, 2 - always
.
   declare $rules=1;
   $ch->add('$rules', <<'.');
# Display the information about the rules:
# 0 - nothing, 1 - significant failures, 2 - summary and all failed preconditions, 3 - executed rule executed
.
   declare $scheduler=0;
   $ch->add('$scheduler', <<'.');
# Reveal the internals of the rule scheduler:
# 0 - nothing, 1 - summary and statistics, 2 - initial rule selection,
# 3 - shortest path search (overwhelming amount of data)
.
   declare $cpp=0;
   $ch->add('$cpp', <<'.');
# Tell about the actions of the perl/C++ interface generator:
# 0 - nothing, 1 - compiler calls and source file updates, 2 - source code generated
.
   declare $files=1;
   $ch->add('$files', <<'.');
# Notify about nontrivial actions during data file processing
.
   declare $external=0;
$ch->add('$external', <<'.');
# Notify about external programs starting in the background
# (not to be mixed up with credits!)
.

   package Polymake::User;

   declare @start_applications;
   $ch->add('@start_applications', <<'.');
# Applications to be loaded at the beginning of each interactive or batch session
.
   declare $default_application;
   $ch->add('$default_application', <<'.');
# Application to start with as the current one
.
   declare @extensions;
   $ch->add('@extensions', <<'.', Core::Customize::State::accumulating);
# A list of directories containing imported and/or locally created extensions
.
   declare %disabled_extensions;
   $ch->add('%disabled_extensions', <<'.', Core::Customize::State::config | Core::Customize::State::hidden | Core::Customize::State::noexport);
# Extensions which could not be configured for given architecture
.
   declare @lookup_scripts;
   $ch->add('@lookup_scripts', <<'.', Core::Customize::State::accumulating);
# A list of directories where to look for scripts
.
   declare $history_size=200;
   $ch->add('$history_size', <<'.');
# Maximal number of commands stored in the interactive shell's history.
# If set to undef, history list grows unlimited.
.
   declare $history_editor=$ENV{VISUAL} || $ENV{EDITOR} || "vi";
   $ch->add('$history_editor', <<'.');
# Editor for the ''history'' command.
# Must be a complete shell command. If the temporary file name is expected somewhere in the middle
# of the arguments, please use the placeholder %f.
.
   declare $help_key="_k1";
$ch->add('$help_key', <<'.');
# Key to press for interactive help in the shell.  Defaults to F1.
.
   declare $help_delimit=1;
$ch->add('$help_delimit', <<'.');
# Add delimiters for better readability in help text.
.

   $ch->cleanup;

   # rescue the old-fashioned lookup_applications list
   if (defined (my $lookup_apps=*lookup_applications{ARRAY})) {
      push @extensions, @$lookup_apps;
      $ch->set('@extensions');
   }

   # treat relative paths as starting at $HOME
   s{^(?:~/|(?!/))}{$ENV{HOME}/} for @extensions, @lookup_scripts;

   if (!@start_applications) {
      $ch->set('@start_applications', map { m{apps/([^/]+)/rules} } glob("$InstallTop/apps/*/rules/main.rules"));
   }
   if (!$default_application) {
      $ch->set('$default_application', string_list_index(\@start_applications, "polytope")>=0 ? "polytope" : $start_applications[0]);
   }
};

1

# Local Variables:
# cperl-indent-level:3
# indent-tabs-mode:nil
# End:
