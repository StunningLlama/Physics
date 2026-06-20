package electrodynamics;

import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;

import com.codedisaster.steamworks.SteamAPI;
import com.codedisaster.steamworks.SteamException;
import com.codedisaster.steamworks.SteamFriends;
import com.codedisaster.steamworks.SteamFriendsCallback;
import com.codedisaster.steamworks.SteamLibraryLoader;
import com.codedisaster.steamworks.SteamLibraryLoaderGdx;
import com.codedisaster.steamworks.SteamNativeHandle;
import com.codedisaster.steamworks.SteamPublishedFileID;
import com.codedisaster.steamworks.SteamRemoteStorage.PublishedFileVisibility;
import com.codedisaster.steamworks.SteamResult;
import com.codedisaster.steamworks.SteamScreenshots;
import com.codedisaster.steamworks.SteamScreenshotsCallback;
import com.codedisaster.steamworks.SteamUGC;
import com.codedisaster.steamworks.SteamUGC.ItemInstallInfo;
import com.codedisaster.steamworks.SteamUGCCallback;
import com.codedisaster.steamworks.SteamUGCDetails;
import com.codedisaster.steamworks.SteamUGCQuery;
import com.codedisaster.steamworks.SteamUGCUpdateHandle;
import com.codedisaster.steamworks.SteamUserStats;
import com.codedisaster.steamworks.SteamUserStatsCallback;
import com.codedisaster.steamworks.SteamUtils;
import com.codedisaster.steamworks.SteamUtilsCallback;

import electrodynamics.gui.SteamUploadUI;

public class Steam {
	//Cloud
	//Workshop
	//Achievements
	//Screenshot
	public static SteamUGC UGC;
	public static SteamScreenshots Screenshots;
	public static SteamUtils Utils;
	public static SteamUserStats UserStats;
	public static SteamFriends Friends;
	public static Simulation e;
	
	public static String ws_title = "";
	public static String ws_description = "";
	private static SteamUploadUI uploadui;
	
	public static void initialize() {
		if (!BuildFlags.steam_enabled)
			return;
		
		try {
		    SteamLibraryLoader loader = new SteamLibraryLoaderGdx();
		    //SteamLibraryLoader loader = new SteamLibraryLoaderLwjgl3();
		    
		    if (!SteamAPI.loadLibraries(loader)) {
		    	System.out.println("SteamAPI.loadLibraries failed");
			    return;
		    }

		    if (!BuildFlags.steam_debugging) {
		    	if (SteamAPI.restartAppIfNecessary(BuildFlags.steam_app_id)) {
		    		System.out.println("Restarting through steam");
		    		System.exit(0);
		    	}
		    }

		    if (!SteamAPI.init()) {
		    	System.out.println("Steam API did not initialize correctly.");
			    return;
		    }

		    UGC = new SteamUGC(new SteamUGCCallback () {
				@Override
				public void onCreateItem(SteamPublishedFileID publishedFileID, boolean needsToAcceptWLA, SteamResult result) {
					if (result == SteamResult.OK) {
						if (needsToAcceptWLA) {
							//TODO
						}

						File folder = SemiSim.getUserFile("_tmp");
						if (!folder.exists()) {
							folder.mkdir();
						}

						File newfile = SemiSim.getUserFile("_tmp/workshop_item.semisim");
						e.savemanager.writeFile(newfile);

			    		File outputimgfile = SemiSim.getUserFile("_tmp/image.png");
				    	try {
				    		if (!outputimgfile.exists())
				    			outputimgfile.createNewFile();
				    		ImageIO.write(e.renderer.img_front, "png", outputimgfile);
				    	} catch (IOException e) {
				    		e.printStackTrace();
				    	}

						SteamUGCUpdateHandle handle = UGC.startItemUpdate(Utils.getAppID(), publishedFileID);
						UGC.setItemVisibility(handle, PublishedFileVisibility.Public);
						UGC.setItemTitle(handle, ws_title);
						UGC.setItemContent(handle, folder.getAbsolutePath());
						UGC.setItemPreview(handle, outputimgfile.getAbsolutePath());
						UGC.setItemDescription(handle, ws_description);
						UGC.submitItemUpdate(handle, "");
					} else {
						SemiSim.displayErrorMessage(new Exception("Steam was not able to create the workshop item."));
					}
				}
				
				@Override
				public void onSubmitItemUpdate(SteamPublishedFileID publishedFileID,
												boolean needsToAcceptWLA, SteamResult result) {
					SemiSim.getUserFile("_tmp/workshop_item.semisim").delete();
		    		SemiSim.getUserFile("_tmp/image.png").delete();
					SemiSim.getUserFile("_tmp").delete();
					if (result == SteamResult.OK) {
						JOptionPane.showMessageDialog(e.opts, "Workshop item successfully uploaded! Redirecting to workshop page...");

						try {
							java.awt.Desktop.getDesktop().browse(new URI("steam://url/CommunityFilePage/" + SteamNativeHandle.getNativeHandle(publishedFileID)));
						} catch (IOException | URISyntaxException ex) {
							ex.printStackTrace();
						}
					}
					else {
						SemiSim.displayErrorMessage(new Exception("Steam was not able to update the workshop item."));
					}
				}

				@Override
				public void onDownloadItemResult(int appID, SteamPublishedFileID publishedFileID, SteamResult result) {
					System.out.println(publishedFileID);
				}
				
				@Override
				public void onUGCQueryCompleted(SteamUGCQuery query, int numResultsReturned, int totalMatchingResults,
												 boolean isCachedData, SteamResult result) {
					
					List<SteamUGCDetails> list = new ArrayList<SteamUGCDetails>();
					List<String> names = new ArrayList<String>();
					for (int i = 0; i < numResultsReturned; i++)
					{
						SteamUGCDetails details = new SteamUGCDetails();
						UGC.getQueryUGCResult(query, i, details);
						list.add(details);
						names.add(details.getTitle());
					}
					
					UGC.releaseQueryUserUGCRequest(query);
					

					JList<Object> tmplist = new JList<>(names.toArray());

					tmplist.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
					tmplist.setVisibleRowCount(5);
					tmplist.setSelectedValue(Preset.DEFAULT, true);

					JScrollPane scrollPane = new JScrollPane(tmplist);
			        scrollPane.setPreferredSize(new Dimension(300, 250));

					int result2 = JOptionPane.showConfirmDialog(null, scrollPane, "Load workshop item", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);

					if (result2 == JOptionPane.OK_OPTION) {
						int selected = tmplist.getSelectedIndex();

						if (selected != -1) {
							SteamPublishedFileID id = list.get(selected).getPublishedFileID();
							ItemInstallInfo info = new ItemInstallInfo();
							UGC.getItemInstallInfo(id, info);
							File file = Paths.get(info.getFolder(), "workshop_item.semisim").toFile();
							SwingUtilities.invokeLater(() -> {
								e.savemanager.readfile(file);
							});
						}
					}
				}
			});
		    
			Screenshots = new SteamScreenshots(new SteamScreenshotsCallback() {
				@Override
				public void onScreenshotRequested() {
					e.controls.takeScreenshot();
				}
			});
			
			Screenshots.hookScreenshots(true);
			
			Utils = new SteamUtils(new SteamUtilsCallback() {
			});

			UserStats = new SteamUserStats(new SteamUserStatsCallback() {
			});

			Friends = new SteamFriends(new SteamFriendsCallback() {
			});
			
			uploadui = new SteamUploadUI();
			uploadui.setVisible(false);
			
			System.out.println("Overlay " + Utils.isOverlayEnabled());
			System.out.println(Friends.getPersonaName());
			System.out.println("Steam " + SteamAPI.isSteamRunning());
		} catch (SteamException e) {
			e.printStackTrace();
		}
	}

	public static void loadUGC() {
		if (SteamAPI.isSteamRunning()) {
			SteamPublishedFileID[] ids = new SteamPublishedFileID[UGC.getNumSubscribedItems(false)];
			UGC.getSubscribedItems(ids, false);

			Collection<SteamPublishedFileID> list = Arrays.asList(ids);
			SteamUGCQuery query = UGC.createQueryUGCDetailsRequest(list);
			UGC.sendQueryUGCRequest(query);
		}
	}

	public static void setAchievement(String ID) {
		if (SteamAPI.isSteamRunning()) {
			if (!UserStats.isAchieved(ID, false)) {
				UserStats.setAchievement(ID);
				e.renderer.achievement_name = ID;
				e.renderer.achievement_timer = 180;
			} else {
				//System.out.println("Already achieved");
			}
		}
	}

	public static void createWorkshopItem() {
		if (SteamAPI.isSteamRunning()) {
			uploadui.desc.setText(e.description);
			uploadui.text_title.setText("(Simulation title)");
			uploadui.setVisible(true);
			//System.out.println("Here");
			/*synchronized(uploadui) {
				try {
					uploadui.wait(Long.MAX_VALUE);
				} catch (InterruptedException e1) {}
			}
			if (uploadui.result == 1) {
				ws_title = uploadui.text_title.getText();
				ws_description = uploadui.desc.getText();
				UGC.createItem(Utils.getAppID(), WorkshopFileType.Community);
			}*/
		}
	}

	public static void addSteamScreenshot(String path, int width, int height) {
		if (SteamAPI.isSteamRunning()) {
			Steam.Screenshots.addScreenshotToLibrary(path, "", width, height);
		}
	}
}
