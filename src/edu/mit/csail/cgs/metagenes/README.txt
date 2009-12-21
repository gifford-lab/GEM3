
Profile class hierarchy: 

            ---> ImmutableProfile ---> PointProfile
          /
Profile  /
         \
          \
            ---> MutableProfile ---> MetaProfile
            
PointProfile are the two primary implementations of Profile.  PointProfile 
is simple, and associates a (String) "name" and a Point object with each array 
of values.

MetaProfile is a Profile built out of underlying Profile objects, and (has the ability
to) average them together.  MetaProfile also listens to its underlying Profiles, so 
if any of them change, it can change too.  (This means that we could, conceivably,
set up a hierarchy of averaged profiles.  This will be useful in the future, 
I promise.) 

ImmutableProfile and MutableProfile differ only in whether they save, or throw 
away, any ProfileListeners -- since ImmutableProfiles (should) never change, they 
don't need to handle the event/messaging listeners.  MutableProfiles, on the other 
hand, need to correctly save and notify listeners upon changes.  Subclasses of 
MutableProfile can call the simple dispatchChange() method to notify listeners of 
any value change in the profile.  