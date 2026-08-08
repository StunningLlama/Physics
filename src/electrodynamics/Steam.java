// Copyright (c) Brandon Li 2026
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
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

import electrodynamics.gui.SteamDownloadUI;
import electrodynamics.gui.SteamUploadUI;

public class Steam {
	public static SteamUGC UGC;
	public static SteamScreenshots Screenshots;
	public static SteamUtils Utils;
	public static SteamUserStats UserStats;
	public static SteamFriends Friends;
	public static Simulation e;
	
	public static String ws_title = "";
	public static String ws_description = "";
	public static SteamUploadUI uploadui;
	public static SteamDownloadUI downloadui;
	
	private static HashMap<String, String> ach_names = new HashMap<String, String>();
	static
	{
		ach_names.put("NEW_MATERIAL", "Material master");
		ach_names.put("WORKSHOP", "Content creator");
		ach_names.put("CRASH", "Party crasher");
		ach_names.put("MANUAL", "Manual labor");
		ach_names.put("LAG", "Halting problem");
	}
	
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
		    
		    Runtime.getRuntime().addShutdownHook(new Thread() {
			    public void run() { shutdown(); }
			});

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
						
						File outputimgfile = SemiSim.getUserFile("_tmp/image.png");
				    	try {
				    		if (!outputimgfile.exists())
				    			outputimgfile.createNewFile();
				    		ImageIO.write(e.renderer.img_front, "png", outputimgfile);
				    	} catch (IOException e) {
				    		e.printStackTrace();
				    	}
				    	
						File newfile = SemiSim.getUserFile("_tmp/workshop_item.semisim");
						e.savemanager.writeFile(newfile, () -> {
							SteamUGCUpdateHandle handle = UGC.startItemUpdate(Utils.getAppID(), publishedFileID);
							UGC.setItemVisibility(handle, PublishedFileVisibility.Public);
							UGC.setItemTitle(handle, ws_title);
							UGC.setItemContent(handle, folder.getAbsolutePath());
							UGC.setItemPreview(handle, outputimgfile.getAbsolutePath());
							UGC.setItemDescription(handle, ws_description);
							UGC.submitItemUpdate(handle, "");
							JOptionPane.showMessageDialog(e.opts, "Upload is starting. Please wait for upload to finish.");
						});
					} else {
						SemiSim.displayWarningMessage("Error!", "Steam was not able to create the workshop item.");
					}
				}
				
				@Override
				public void onSubmitItemUpdate(SteamPublishedFileID publishedFileID,
												boolean needsToAcceptWLA, SteamResult result) {
					SemiSim.getUserFile("_tmp/workshop_item.semisim").delete();
		    		SemiSim.getUserFile("_tmp/image.png").delete();
					SemiSim.getUserFile("_tmp").delete();
					if (result == SteamResult.OK) {
						Steam.setAchievement("WORKSHOP");
						JOptionPane.showMessageDialog(e.opts, "Workshop item successfully uploaded! Redirecting to workshop page...");

						try {
							java.awt.Desktop.getDesktop().browse(new URI("steam://url/CommunityFilePage/" + SteamNativeHandle.getNativeHandle(publishedFileID)));
						} catch (IOException | URISyntaxException ex) {
							ex.printStackTrace();
						}
					}
					else {
						SemiSim.displayWarningMessage("Error!", "Steam was not able to update the workshop item.");
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
					
					downloadui.list.setModel(new DefaultComboBoxModel<>(names.toArray(new String[] {})));
					downloadui.details = list;
					downloadui.setVisible(true);
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

			downloadui = new SteamDownloadUI();
			downloadui.setVisible(false);
			
			System.out.println("Overlay " + Utils.isOverlayEnabled());
			System.out.println(Friends.getPersonaName());
			System.out.println("Steam " + SteamAPI.isSteamRunning());
		} catch (SteamException e) {
			e.printStackTrace();
		}
	}
	
	public static void loadUGCs() {
		if (SteamAPI.isSteamRunning()) {
			SteamPublishedFileID[] ids = new SteamPublishedFileID[UGC.getNumSubscribedItems(false)];
			UGC.getSubscribedItems(ids, false);

			Collection<SteamPublishedFileID> list = Arrays.asList(ids);
			SteamUGCQuery query = UGC.createQueryUGCDetailsRequest(list);
			UGC.sendQueryUGCRequest(query);
		}
	}

	public static void openUGC(SteamPublishedFileID id) {
		if (SteamAPI.isSteamRunning()) {
			ItemInstallInfo info = new ItemInstallInfo();
			UGC.getItemInstallInfo(id, info);
			File file = Paths.get(info.getFolder(), "workshop_item.semisim").toFile();
			new Thread(() -> {
				e.savemanager.readfile(file);
			}).start();
		}
	}


	public static BufferedImage loadThumbnail(SteamPublishedFileID id) {
		if (SteamAPI.isSteamRunning()) {
			ItemInstallInfo info = new ItemInstallInfo();
			UGC.getItemInstallInfo(id, info);
			File file = Paths.get(info.getFolder(), "image.png").toFile();
			try {
				return ImageIO.read(file);
			} catch (IOException e) {}
		}

		return null;
	}

	public static void setAchievement(String ID) {
		if (SteamAPI.isSteamRunning()) {
			if (!UserStats.isAchieved(ID, false)) {
				UserStats.setAchievement(ID);
				e.renderer.achievement_name = ach_names.get(ID);
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

	public static void shutdown() {
		if (SteamAPI.isSteamRunning()) {
			UGC.dispose();
			Screenshots.dispose();
			Utils.dispose();
			UserStats.dispose();
			Friends.dispose();
			SteamAPI.shutdown();
		}
	}
}
